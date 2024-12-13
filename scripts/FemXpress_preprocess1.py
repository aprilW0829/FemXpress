#!usr/bin/python3
from Bio import SeqIO
import subprocess
import pysam
import os
import sys
import vcf
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
import fileinput
import math
import resource
import seaborn as sns
from optparse import OptionParser
import argparse
import shutil
import codecs
import chardet
import shutil

parser = argparse.ArgumentParser(description='preprocess of FemXpress,\nExample:python FemXpress_1.py -b possorted_genome_bam.bam -g genome.fa -e meta.txt -r rmsk.txt')
parser.add_argument("-b","--bam", type=str,help="input bam file from CellRanger output")
parser.add_argument("-g","--genome",type=str,help="input genome")
parser.add_argument("-e","--meta",type=str,help="input meta information of barcodes")
parser.add_argument("-r","--rmsk",type=str,help="rmsk file")
parser.add_argument("-n","--sample",type=str,help="sample name")
parser.add_argument("-c","--minimum_counts", type=int,help="minimum counts supporting second largest base, default: 2")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=False, help="print verbose output")
args=parser.parse_args()
#print(args.bam)
#print(args.genome)
#print(args.meta)
#print(args.rmsk)

## imput_file
# sys.argv[1]:input_bam (unfiltered bam from cellranger output)
# sys.argv[2]:ref_genome_file (ref.fa)
# sys.argv[3]:meta_data for scRNA-seq(barcodes after celltype annotation)
# sys.argv[4]:rmsk.txt (rmsk txt downloaded from UCSC)
# sample name:first str of bam
## output_file
# sys.argv[5]:final_filtered_snp.txt(after snp and manual filtering based on sequence)

# define locations of used tools
current_dir = os.path.dirname(os.path.realpath(__file__))
samtools = os.path.join(current_dir, '..', 'apps', 'samtools')
bam_readcount = os.path.join(current_dir, '..', 'apps', 'bam-readcount')
freebayes = os.path.join(current_dir, '..', 'apps', 'freebayes')

# set to incerease txt number at the same time
soft_limit, hard_limit = resource.getrlimit(resource.RLIMIT_NOFILE)
print("[FemXpress preprocess] file number's limit:", soft_limit, hard_limit)
new_limit = 63655
resource.setrlimit(resource.RLIMIT_NOFILE, (new_limit, hard_limit))
soft_limit, hard_limit = resource.getrlimit(resource.RLIMIT_NOFILE)
#print("new txt number's limit:", soft_limit, hard_limit)

pd.set_option('display.max_rows', 20)
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_colwidth', 200)

# make tmp dir and result dir
tmp_dir='./FemXpress/'+args.sample+'/tmp'
result_dir='./FemXpress/'+args.sample+'/result'
inter_base_dir='./FemXpress/'+args.sample+'/base_count'


################################################# preprocess: cellranger and get passorted_chrX.bam(unfiltered),but nedd to do celltypes' annotation by yourself
# define function that opens files
def get_encoding(file):
    with open(file,'rb') as f:
        return chardet.detect(f.read())['encoding']

# define function to create a new directory
def create_three_inter_dir(dir_paths):
#dir_paths=[tmp_dir,result_dir,inter_base_dir]
    for dir_path in dir_paths:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            print('[FemXpress preprocess] create tmp dir!')
        else:
            print('[FemXpress preprocess] tmp dir already exist!')

# define function thqt can get all files under the indicated dir
#all_d_snp_loci_dic={}
def all_path(dirname,defined_txt_str):
    sgtxt_list = []
    for file_name_list in os.walk(dirname,defined_txt_str):
        for filename in file_name_list:
            for i in filename:
                if i.endswith(defined_txt_str):
                    sgtxt_list.append(i)
    return sgtxt_list

################################################# step1: get each barcode's bam from 10×Genomics
# define function that can get index for chrX.bam
def get_index_for_bam(bam):
    subprocess.run([samtools, "index", bam])
    print('[FemXpress preprocess] create index for all.bam!')

# define function that can get snp for chrX.bam using freebayes using the following parameters
def get_chrX_bam(input_bam,output_bam):
    global samtools
    global tmp_dir
    output_bam=tmp_dir+'/output_chrX.bam'
    chromosome='chrX'
    command2 = [samtools, "view", "-b", "-h", input_bam, chromosome]
    with open(output_bam, "w") as output:
        process = subprocess.Popen(command2, stdout=output)
        process.wait()
    print('[FemXpress preprpcess] get chrX.bam!')
    subprocess.run([samtools, "index", output_bam])
    print('[FemXpress preprocess] create index for chrX.bam!')

# define function that can read metadata (have filtered using Seurat, scanpy, et al., or not)
def get_barcodes(f):
    #f_meta=open(args.meta,'r')
    global tmp_dir
    trimmed_barcode_list = []
    f_meta=open(f,'r')
    fo_tmp1=open(tmp_dir+'/all_barcodes.txt','w')
    start_pos = None
    end_pos = None
    barcode_list=[]
    trimmed_barcode_list=[]
    for i in f_meta:
        i=i.split('\t')
        barcode_list.append(i[0].strip('\n'))
    for barcode in barcode_list:
        ###search for first location (up than ACGT four bases)
        for i in range(len(barcode)-16):
            if set(barcode[i:i+16]).issubset({'A','C','G','T'}):
                start_pos=i
                break
        if start_pos is not None:
            end_pos = -1
        trimmed_barcode = barcode[start_pos:]
        fo_tmp1.write(trimmed_barcode+'\n')
        trimmed_barcode_list.append(trimmed_barcode)
    print('[FemXpress preprocess] add all barcodes to meta info,done!')
    #print(len(trimmed_barcode_list))
    return trimmed_barcode_list

# define function that can get chrX.bam of indicated barcodes
def get_filtered_chrX_bam(input_bam, output_bam, barcode_list):
    global trimmed_barcode_list
    global tmp_dir
    #input_bam=tmp_dir+'/output_chrX.bam'
    #output_bam=tmp_dir+'/filtered_chrX.bam'
    #print(len(barcode_list))
    bam_in=pysam.AlignmentFile(input_bam, "rb")
    bam_out=pysam.AlignmentFile(output_bam, "wb", template=bam_in)
    #chrX;6054451    0;1;0;0 C       1
    for read in bam_in.fetch(until_eof=True):
        if read.has_tag('CB'):
            cb_tags=dict(read.tags)['CB']
            if cb_tags in barcode_list:
                bam_out.write(read)
    print('[FemXpress preprocess] get filtered bam done!')

############################# get final vcf for download analysis(SNP and manual filtering based on sequence)
# define function that can get vcf of filtered_chrX.bam
def get_vcf_of_filtered_chrX_bam(ref_file, bam_file, output_vcf):
    #global freebayes
    global tmp_dir
    #ref_file=args.genome
    #bam_file=tmp_dir+'/filtered_chrX.bam'
    #output_vcf=tmp_dir+'/filtered_freebayes.vcf'
    print('[FemXpress preprocess] call SNP using freebyes,start!')
    #cmd2 = f"{freebayes} -f {ref_file} {bam_file} -C 1 > {output_vcf}"
    #output=subprocess.check_output(cmd2,shell=True)
    #subprocess.run([freebayes, "-f" ,ref_file ,bam_file ,"-C", "1", ">", output_vcf])
    cmd = [freebayes, "-f", ref_file, bam_file, "-C 1"]
    with open(output_vcf, 'w') as f_out:
        result = subprocess.run(cmd, check=True, stdout=f_out, text=True)
    print('[FemXpress preprocess] get vcf for filtered chrX bam done!')

# define function that can get filtered SNP from the whole vcf
def get_filtered_snp_from_vcf(input_vcf):
    global tmp_dir
    #input_vcf=tmp_dir+'/filtered_freebayes.vcf'
    vcf_reader = vcf.Reader(filename=input_vcf)
    output_vcf = vcf.Writer(open(tmp_dir+"/filtered_freebayes_snps.vcf", "w"), vcf_reader)
    sample_name = vcf_reader.samples[0]
    for record in vcf_reader:
        genotype = record.genotype(sample_name)['GT']
        if record.is_snp:
            if genotype == '1/1' or genotype == '0/0' or genotype == './.':
                pass
            else:
                output_vcf.write_record(record)
    output_vcf.close()

# define function that can get bed of filtered SNP
def get_bed_from_filtered_snp(input_vcf, bed_file):
    global tmp_dir
    #input_vcf=tmp_dir+'/filtered_freebayes_snps.vcf'
    #bed_file = open(tmp_dir+'/filtered_freebayes_snps.bed', 'w')
    vcf_reader = vcf.Reader(filename=input_vcf)
    for record in vcf_reader:
        chrom = record.CHROM
        start = record.POS
        end = start + 1
        bed_file.write(f"{chrom}\t{start}\t{end}\n")
    bed_file.close()

# define function that can filter SNP based on artificial conditions
# define function that can get chrX of indicated genome(such as: mm10, hg38, et al.)
def get_chrX_fa(genome_file):
    seq1='ACGTACGTACGTACGT'
    #genome_file = args.genome
    genome = SeqIO.parse(genome_file, "fasta")
    for sequence in genome:
        if sequence.id == 'chrX':
            seq1 = sequence.seq
    #print(seq1)
    print('[FemXpress preprocess] get fasta of chrX done!')
    return seq1

# define function that can onstruct simple repeat loci dic
def get_simple_repeat_snp(rmsk_file):
    simple_repeat_snp_dic = {}
    f=open(rmsk_file,'r')
    for a in f:
        a=a.split('\t')
        if a[11] =='Simple_repeat':
            chr_seq=a[5]
            for b in range(int(a[6]),int(a[7].strip('\n'))+1):
                simple_loci=chr_seq+';'+str(b)
                simple_repeat_snp_dic[simple_loci]=0
    print('[FemXpress preprocess] get simple repeat sequence of mm10"s chrX of mm10 done!')
    return simple_repeat_snp_dic


# get up and down 20bp bed location and sequence of the SNP (filter simple_repeat_sequence :yes and or)
#def get_up_down20bp_bed_sequence_dic(f,genome1,simple_repeat_snp_dic):
# imput: f1=open(tmp_dir+'/filtered_freebayes_snps.bed','r')
#                           simple_repeat_snp_dic
def get_filtered_snp_2round(input_file, dic, seq1):
    snp_updown20bp_bed_sequence_dic={}
    f = open(input_file, 'r')
    for i in f.readlines():
        i=i.split('\t')
        chr_seq=i[0]
        loci=int(i[1])
        loci1=loci+1
        snp_loci=chr_seq+';'+str(loci)
        start=loci-10
        end=loci+11
        bed_info=chr_seq+'\t'+str(start)+'\t'+str(end)
        snp_sequence=seq1[start:end].upper()
        snp_seq_info=chr_seq+'\t'+str(loci)+'\t'+str(loci1)+'\t'+snp_sequence
        #snp_updown20bp_bed_sequence_ml_dic[snp_loci]=snp_sequence
        if snp_loci in dic.keys():
            #print(snp_loci)
            pass
        else:
            if snp_seq_info in snp_updown20bp_bed_sequence_dic.keys():
                pass
            else:
                snp_updown20bp_bed_sequence_dic[snp_seq_info]=0
    return snp_updown20bp_bed_sequence_dic


#                           snp_updown20bp_bed_sequence_dic                       
def get_filtered_snp_3round(dic):
    fo_dic = {}
    fo1_dic = {}
    fo2_dic = {}
    fo3_dic = {}
    fo3 = open(tmp_dir + '/filtered_freebayes_snps_manual_info.txt', 'w')

    ### 1) remove sequence that includes only two bases.
    all_line_dic={}
    all_base_dic={}
    all_base_list=[]
    for i in dic.keys():
        x=i
        i=i.split('\t')
        sequence=str(i[3].strip('\n'))
        for j in range(len(sequence)):
            a=sequence[j].upper()
            if a in all_base_dic.keys():
                all_base_dic[a]+=1
            else:
                all_base_dic[a]=0
                all_base_list.append(a)
        if len(all_base_list)<=2:
            pass
        else:
            all_line_dic[x.strip('\n')]=0
            fo_dic[x]=0
        all_base_dic={}
        all_base_list=[]
        all_line_dic={}
    ### 2) remove sequence,example:AAAAAAAAA,TTTTTTTT et al.
    for i in fo_dic.keys():
        x=i
        i=i.split('\t')
        #sequence=i[18].upper()
        sequence=i[3].strip('\n')
        sequence=str(sequence).upper()
        #print(type(sequence))
        pattern1=re.compile(r'(A{5,21})')
        pattern2=re.compile(r'(C{5,21})')
        pattern3=re.compile(r'(G{5,21})')
        pattern4=re.compile(r'(T{5,21})')
        r1=re.search(pattern1,sequence)
        r2=re.search(pattern2,sequence)
        r3=re.search(pattern3,sequence)
        r4=re.search(pattern4,sequence)
        if r1:
            pass
        elif r2:
            pass
        elif r3:
            pass
        elif r4:
            pass
        else:
            fo2_dic[x]=0
    ### 3) remove sequnece including seq that longer than 10bp
    A_num=0
    C_num=0
    G_num=0
    T_num=0
    seq_list =[]
    for i in fo2_dic.keys():
        x=i
        i=i.split('\t')
        sequence=str(i[3].strip('\n')).upper()
        seq_list = list(sequence)
        A_num = seq_list.count('A')
        C_num = seq_list.count('C')
        G_num = seq_list.count('G')
        T_num = seq_list.count('T')
        if A_num >= 10:
            pass
        elif C_num >= 10:
            pass
        elif G_num >= 10:
            pass
        elif T_num >= 10:
            pass
        else:
            #print(x)
            fo3_dic[x]=0
            #fo3.write(x)
    for a in fo3_dic.keys():
        #print(a.strip('\n'))
        fo3.write(str(a)+'\n')
    return fo3

#src_file = tmp_dir+'/filtered_freebayes_snps.bed'
#dst_file = tmp_dir+'/filtered_freebayes_snps_manual_info.txt'
#shutil.copy(src_file, dst_file)

############################################ single barcode's bam split by barcode from filtered_chrX.bam,filtering snp by judge base counts in all barcodes' ratio,plotting unfiltered snps' distribution and filtered snps' distribution,get final snp-barcode matrix
# define function that can get each barcode's filtered chrX.bam
def get_split_bam(input_bam):
    global inter_base_dir
    #input_bam=tmp_dir+'/filtered_chrX.bam'
    output_dir = inter_base_dir
    with pysam.AlignmentFile(input_bam, "rb",ignore_truncation=True) as bam_file:
        barcode_files = {}
        for read in bam_file:
            barcode = read.get_tag("CB")
            if barcode not in barcode_files:
                filename = f"{output_dir}/{barcode}.bam"
                #print(barcode)
                outfile = pysam.AlignmentFile(filename, "wb", template=bam_file)
                barcode_files[barcode] = outfile
            barcode_files[barcode].write(read)
    for outfile in barcode_files.values():
        outfile.close()


############################################# step2: get each snp's ACGT count for each barcode
def get_each_snp_ACGT_count_for_each_barcode(bam_files, bed_file, genome):
    
    global inter_base_dir
    global tmp_dir
    
    base_count_files = []
    for bam_file in bam_files:
        bam_file1 = inter_base_dir + '/' + bam_file
        barcode = bam_file[:-4]
        #print(barcode)
        base_count_file1 = inter_base_dir + '/' + barcode + '_base_count1.txt'
        # chrX    6407607 T       1       =:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  A:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00     G:1:255.00:25.00:255.00:0:1:0.16:0.00:25.00:0:0.00:150.00:0.08  T:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00

        #command = f"{bam_readcount} -w 0 -f {genome} {bam_file1} -l {bed_file} | awk 'BEGIN{{FS=OFS=\"\t\"}}{{split($6,A,\":\");split($7,C,\":\");split($8,G,\":\");split($9,T,\":\");print $1\";\"$2,A[2]\";\"C[2]\";\"G[2]\";\"T[2],$3,$4}}' > {base_count_file}"
        command = "{} -w 0 -f {} {} -l {} > {}".format(bam_readcount, genome, bam_file1, bed_file, base_count_file1)
        os.system(command)
        #subprocess.run(command, shell=True)

        if os.path.getsize(base_count_file1) == 0:
            pass
        else:
            base_count_file = open(inter_base_dir + '/' + barcode + '_base_count.txt', 'w')
            base_count_file.write('loci'+'\t'+'A;C;G;T'+'\t'+'ref'+'\t'+'total_count'+'\n')
            base_count_file_pre = open(inter_base_dir + '/' + barcode + '_base_count1.txt', 'r')
            for i in base_count_file_pre:
                i = i.split('\t')
                chr_seq = i[0]
                snp_loci = i[1]
                A_count = i[5].split(':')[1]
                C_count = i[6].split(':')[1]
                G_count = i[7].split(':')[1]
                T_count = i[8].split(':')[1]
            
                #print(chr_seq+';'+snp_loci+'\t'+A_count+';'+C_count+';'+G_count+';'+T_count+'\t'+i[2].upper()+'\t'+i[3])
                base_count_file.write(chr_seq+';'+snp_loci+'\t'+A_count+';'+C_count+';'+G_count+';'+T_count+'\t'+i[2].upper()+'\t'+i[3]+'\n')

                base_count_files.append(base_count_file)
        # print('get '+barcode+' bam base count done!')
    #print(len(base_count_files))
    print('[FemXpress preprocess] get all barcode bam\'s base count done!')
    return base_count_files


# define function that run pipeline
def runpipeline(a):
    global tmp_dir
    global inter_base_dir
    global result_dir
    global samtools
    global trimmed_barcode_list
    global sgtxt_list

    ############################# preprocess: cellranger and get passorted_chrX.bam(unfiltered),but nedd to do celltypes' annotation by yourself
    # create temporary directories
    dir_paths = [tmp_dir, result_dir, inter_base_dir]
    create_three_inter_dir(dir_paths)
    
    # index for bam
    bam =args.bam
    get_index_for_bam(bam)
    
    # get chrX.bam of possorted_genome_bam.bam
    bam =args.bam
    output_bam = tmp_dir + '/output_chrX.bam'
    get_chrX_bam(bam, output_bam)

    # read barcodes for bam (can be filtered by using some parameters, or not)
    f =args.meta
    trimmed_barcode_list = get_barcodes(f)

    # get chrX.bam of filtered barcodes (all barcodes or filtered barcodes using some quality control conditions)
    input_bam =tmp_dir + '/output_chrX.bam'
    output_bam =tmp_dir + '/filtered_chrX.bam'
    get_filtered_chrX_bam(input_bam, output_bam, trimmed_barcode_list)


    ############################## step1: get each barcode's bam from 10×Genomics
    # index for filtered_chrX.bam
    bam = tmp_dir + '/filtered_chrX.bam'
    get_index_for_bam(bam)
    
    # get vcf for chrX.bam using freebayes
    ref_file=args.genome
    bam_file=tmp_dir+'/filtered_chrX.bam'
    output_vcf=tmp_dir+'/filtered_freebayes.vcf'
    get_vcf_of_filtered_chrX_bam(ref_file,bam_file,output_vcf)

    # get filtered vcf for chrX.bam (Round1 filtering: SNP,0/1, not SNV et al.)
    input_vcf=tmp_dir+'/filtered_freebayes.vcf'
    get_filtered_snp_from_vcf(input_vcf)

    # get bed format of filtered vcf for chrX.bam of Round1 filtering
    input_vcf = tmp_dir + '/filtered_freebayes_snps.vcf'
    bed_file = open(tmp_dir + '/filtered_freebayes_snps.bed', 'w')
    get_bed_from_filtered_snp(input_vcf, bed_file)
    
    # get sequence of chrX of mm10
    genome_file = args.genome
    seq1 = get_chrX_fa(genome_file)

    # get rmsk locations of mm10
    rmsk_file = args.rmsk
    simple_repeat_snp_dic = get_simple_repeat_snp(rmsk_file)

    # filtered vcf for chrX.bam of Round2 filtering (filtering repeat mask)
    input_file =tmp_dir+'/filtered_freebayes_snps.bed'
    genome = args.genome
    snp_updown20bp_bed_sequence_dic = get_filtered_snp_2round(input_file, simple_repeat_snp_dic, seq1)

    # filtered vcf for chrX.bam of Round3 filtering (filtering using customized standards)
    fo3 = get_filtered_snp_3round(snp_updown20bp_bed_sequence_dic)

    # get ach barcode's filtered chrX.bam
    input_bam = tmp_dir + '/filtered_chrX.bam'
    get_split_bam(input_bam)

    # index for each barcode's filtered chrX.bam
    sgtxt_list = []
    sgtxt_list = all_path(inter_base_dir, 'bam')
    print('[FemXpress preprocess] length of split bams is '+str(len(sgtxt_list)))
    for i in sgtxt_list:
        small_bam_file = inter_base_dir + '/' + i
        subprocess.run([samtools, "index", small_bam_file])
    print('[FemXpress preprocess] create index for all the barcodes, done!')
    

    ############################# step2: get ACGT count of all analysed SNPs for each barcode
    sgtxt_list = []
    sgtxt_list = all_path(inter_base_dir, '.bam')
    bam_files = sgtxt_list
    genome = args.genome
    bed_file = tmp_dir + '/filtered_freebayes_snps_manual_info.txt'
    bed_file2 = tmp_dir+'/filtered_freebayes_snps.bed'
    print('[FemXpress preprocess] length of analysed split bams is '+str(len(bam_files)))
    base_count_files = get_each_snp_ACGT_count_for_each_barcode(bam_files, bed_file, genome)
    print('[FemXpress preprocess] length of ACGT count file for bams is '+str(len(base_count_files)))
    

    print('[FemXpress preprocess] preprocess data, step1,done!')

x=1
runpipeline(x)
