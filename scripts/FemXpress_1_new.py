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
print("txt number's limit:", soft_limit, hard_limit)
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
            print('create tmp dir!')
        else:
            print('tmp dir already exist!')

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
    print('create index for all.bam!')

# define function that can get snp for chrX.bam using freebayes using the following parameters
def get_chrX_bam(input_bam,output_bam):
    global samtools
    global tmp_dir
    output_bam=tmp_dir+'/output_chrX.bam'
    chromosome='chrX'
    command2 = [samtools, "view", "-b", "-h", input_bam, chromosome]
    with open(utput_bam, "w") as output:
        process = subprocess.Popen(command2, stdout=output)
        process.wait()
    print('get chrX.bam!')
    subprocess.run([samtools, "index", output_bam])
    print('create index for chrX.bam!')

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
    print('add all barcodes to meta info,done!')
    print(len(trimmed_barcode_list))
    return trimmed_barcode_list

# define function that can get chrX.bam of indicated barcodes
def get_filtered_chrX_bam(input_bam, output_bam, barcode_list):
    global trimmed_barcode_list
    global tmp_dir
    #input_bam=tmp_dir+'/output_chrX.bam'
    #output_bam=tmp_dir+'/filtered_chrX.bam'
    print(len(barcode_list))
    bam_in=pysam.AlignmentFile(input_bam, "rb")
    bam_out=pysam.AlignmentFile(output_bam, "wb", template=bam_in)
    #chrX;6054451    0;1;0;0 C       1
    for read in bam_in.fetch(until_eof=True):
        if read.has_tag('CB'):
            cb_tags=dict(read.tags)['CB']
            if cb_tags in barcode_list:
                bam_out.write(read)
    print('get filtered bam done!')

############################# get final vcf for download analysis(SNP and manual filtering based on sequence)
# define function that can get vcf of filtered_chrX.bam
def get_vcf_of_filtered_chrX_bam(ref_file, bam_file, output_vcf):
    #global freebayes
    global tmp_dir
    #ref_file=args.genome
    #bam_file=tmp_dir+'/filtered_chrX.bam'
    #output_vcf=tmp_dir+'/filtered_freebayes.vcf'
    print('call SNP using freebyes,start!')
    #cmd2 = f"{freebayes} -f {ref_file} {bam_file} -C 1 > {output_vcf}"
    #output=subprocess.check_output(cmd2,shell=True)
    #subprocess.run([freebayes, "-f" ,ref_file ,bam_file ,"-C", "1", ">", output_vcf])
    cmd = [freebayes, "-f", ref_file, bam_file, "-C 1"]
    with open(output_vcf, 'w') as f_out:
        result = subprocess.run(cmd, check=True, stdout=f_out, text=True)
    print('get vcf for filtered chrX bam done!')

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
    print('get fasta of chrX done!')
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
    print('get simple repeat sequence of mm10"s chrX of mm10 done!')
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
    global bam_readcount
    #sgtxt_list = []
    #all_path(inter_base_dir, '.bam')
    #bam_files = sgtxt_list
    #genome = args.genome
    #bed_file = tmp_dir + '/filtered_freebayes_snps_manual_info.txt'
    #bed_file2=tmp_dir+'/filtered_freebayes_snps.bed'
    for bam_file in bam_files:
        bam_file1 = inter_base_dir + '/' + bam_file
        barcode = bam_file[:-4]
        command = f"{bam_readcount} -w 0 -f {genome} {bam_file1} -l {bed_file} | awk 'BEGIN{{FS=OFS=\"\t\"}}{{split($6,A,\":\");split($7,C,\":\");split($8,G,\":\");split($9,T,\":\");print $1\";\"$2,A[2]\";\"C[2]\";\"G[2]\";\"T[2],$3,$4}}' > {inter_base_dir}/{barcode}_base_count.txt"
        subprocess.run(command, shell=True)
        # print('get '+barcode+' bam base count done!')
    print('get all barcode bam\'s base count done!')


############################################ step3: get merged ACGT count of all analysed SNPs for each barcode
# define dunction that can get equal snps for each barcode (using whole snps, can be ,;,;,;,)
def get_merged_snps_ACGT_count_file_for_each_barcode(manual_file):
    data = [] 
    positions = []
    global inter_base_dir
    global tmp_dir
    #with open(tmp_dir+'/filtered_freebayes_snps_manual_info.txt', 'r') as bed_file:
    with open(manual_file,'r') as bed_file:
        for line in bed_file:
            chrom, start, end,up_down_20bp_sequence = line.strip('\n').split('\t')
            #chrom, start, end = line.strip('\n').split('\t')
            a=str(chrom)+';'+str(start)
            data.append(a)
            positions.append((chrom,int(start)))
    data=list(set(data))
    print('number of snp analysed is '+str(len(positions)))
    return data

    columns_names=['loci','A;C;G;T','ref','total_count']
    sgtxt_list=[]
    all_path(inter_base_dir,'-1_base_count.txt')
    count_files=sgtxt_list
    for count_file in count_files:
        #with codecs.open(inter_base_dir+'/'+count_file,encoding='gbk') as file:
        #df = pd.read_csv(file, sep='\t')
        #error_bad_lines = False
        df = pd.read_csv(inter_base_dir+'/'+count_file,names=columns_names,sep='\t')
        fw = open(inter_base_dir+'/'+count_file[:18]+'_merged_base_count.txt','w')
        df = df.drop_duplicates(subset=['loci'])
        search_col = 'loci'
        output_col = 'A;C;G;T'
        for search_value in data:
            if search_value in df[search_col].values:
                row=df.loc[df[search_col]==search_value]
                result1 = np.array2string(row[['loci','A;C;G;T']].values,separator='\t', formatter={'int_kind':lambda x: "%.2f" % x})
                result1 = result1.replace('[','').replace(']','').replace("'","").strip("'")
                if result1.split('\t')[1]=='0;0;0;0':
                    fw.write(result1.split('\t')[0]+ '\t.;.;.;.\n')
                elif result1.split('\t')[1]==' ':
                    fw.write(result1.split('\t')[0]+ '\t.;.;.;.\n')
                else:
                    fw.write(result1+ '\n')
            else:
                fw.write(search_value+'\t.;.;.;.\n')
        #print('get merged_base_count for '+count_file[:18]+' done!')
    print('get merged_base_count for each barcode done!')


############################################# step4: get snp_information
def define_filtered_snp_loci(input_file):
    f=open(input_file,'r')
    snp_loci_list = []
    snp_loci_dic = {}
    # f=open(tmp_dir+'/filtered_freebayes_snps_manual_info.txt','r')
    for i in f:
        i = i.split('\t')
        snp_loci = i[0] + ';' + i[1]
        snp_loci_dic[snp_loci] = 1
        snp_loci_list.append(snp_loci)
    return snp_loci_list, snp_loci_dic

def get_single_snp_coverage(txt_list,snp_loci_dic,snp_loci_list1):
    A_num=0
    C_num=0
    G_num=0
    T_num=0
    ## cell barcode frequency dic
    each_snp_info_dic={}
    ## base count frequency dic
    each_snp_info_dic1={}
    each_snp_info_list=[]
    each_snp_base_list=['A','C','G','T']
    each_snp_base_num=[0,0,0,0]
    each_snp_base_num1=[0,0,0,0]
    snp_loci_4base_barcodes_list_dic = {}

    # chrX;8144585    2;0;0;0 A       2
    snp_loci_4base_barcodes_list_dic = {key: {} for key in snp_loci_list1}
    for a in snp_loci_4base_barcodes_list_dic.keys():
        snp_loci_4base_barcodes_list_dic[a]['A'] = []
        snp_loci_4base_barcodes_list_dic[a]['C'] = []
        snp_loci_4base_barcodes_list_dic[a]['G'] = []
        snp_loci_4base_barcodes_list_dic[a]['T'] = []

    for i in txt_list:
        barcode=i[:18]
        f=open(inter_base_dir+'/'+i,'r')
        for j in f:
            each_snp_base_num=[0,0,0,0]
            each_snp_base_num1=[0,0,0,0]
            each_snp_info_list=[]
            j=j.split('\t')
            snp_loci=j[0]
            if snp_loci in snp_loci_dic.keys():
                if j[1].strip('\n') != '.;.;.;.' :
                    if j[1].strip('\n') != '0;0;0;0':
                        if j[1].strip('\n') != 'A_num;C_num;G_num;T_num':
                            A_num=int(j[1].strip('\n').split(';')[0])
                            C_num=int(j[1].strip('\n').split(';')[1])
                            G_num=int(j[1].strip('\n').split(';')[2])
                            T_num=int(j[1].strip('\n').split(';')[3])
                            each_snp_info_list.extend((A_num,C_num,G_num,T_num))
                            if len(np.nonzero(each_snp_info_list)[0])==1:
                                nonzero_index=[a for a,e in enumerate(each_snp_info_list) if e != 0]
                                nonzero_index_num=nonzero_index[0]
                                #print(nonzero_index_num)
                                nonzero_index_value=each_snp_info_list[nonzero_index_num]
                                #print(str(nonzero_index_num)+'\t'+str(nonzero_index_value))
                                nonzero_base=each_snp_base_list[nonzero_index_num]
                                if snp_loci in each_snp_info_dic.keys() and snp_loci in each_snp_info_dic1.keys():
                                    value=each_snp_info_dic[snp_loci]
                                    value1=each_snp_info_dic1[snp_loci]
                                    value[nonzero_index_num]+=1
                                    value1[nonzero_index_num]+=nonzero_index_value
                                    #print(str(value)+'\t'+str(value1))
                                    each_snp_info_dic[snp_loci]=value
                                    #print(str(each_snp_info_dic[snp_loci]))
                                    each_snp_info_dic1[snp_loci]=value1
                                    #print(str(each_snp_info_dic1[snp_loci]))
                                else:
                                    each_snp_base_num[nonzero_index_num]+=1
                                    each_snp_info_dic[snp_loci]=each_snp_base_num
                                    each_snp_base_num1[nonzero_index_num]+=nonzero_index_value
                                    each_snp_info_dic1[snp_loci]=each_snp_base_num1

                                string='\t'.join('%s' %id for id in each_snp_info_dic[snp_loci])
                                string1='\t'.join('%s' %id for id in each_snp_info_dic1[snp_loci])
                                #print(string1)

                                if nonzero_base=='A':
                                    snp_loci_4base_barcodes_list_dic.setdefault(snp_loci,{})['A'].append(barcode)
                                if nonzero_base=='C':
                                    snp_loci_4base_barcodes_list_dic.setdefault(snp_loci,{})['C'].append(barcode)
                                if nonzero_base=='G':
                                    snp_loci_4base_barcodes_list_dic.setdefault(snp_loci,{})['G'].append(barcode)
                                if nonzero_base=='T':
                                    snp_loci_4base_barcodes_list_dic.setdefault(snp_loci,{})['T'].append(barcode)

    return each_snp_info_dic,each_snp_info_dic1,snp_loci_4base_barcodes_list_dic


################################################# step5: get two main bases's distribution of each snp (unfiltered)
# chrX    4060222 4060223 TTTCTCTGTATAGCCCTGGCT
def get_final_snp_dic(input_file):
    snp_sequence_dic = {}
    #input_file=open(tmp_dir+'/filtered_freebayes_snps_manual_info.txt','r')
    f=open(input_file,'r')
    for i in f:
        i=i.split('\t')
        snp_loci=i[0]+';'+i[1]
        snp_sequence_dic[snp_loci]=i[3].strip('\n')
    #print(len(snp_sequence_dic))
    return snp_sequence_dic

# define function that get reads and counts for unfiltered SNPs
def get_unfiltered_reads(input_file):
    reads_value_list = []
    reads_value_ratio = []
    reads_value_ratio_dic = {}
    df2 = open(input_file,'r')
    for i in df2:
        i=i.split('\t')
        reads=i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4].strip('\n')
        reads_value_list.extend((i[1],i[2],i[3],i[4].strip('\n')))
        reads_value_list = list(map(lambda x: int(x), reads_value_list))
        #print(reads_value_list)
        max_index = reads_value_list.index(max(reads_value_list))
        second_max_index = reads_value_list.index(sorted(reads_value_list)[-2])
        max_value = reads_value_list[max_index]
        second_max_value = reads_value_list[second_max_index]
        #print(str(max_value)+'\t'+str(second_max_value))
        if second_max_value != 0:
            if max_index > second_max_index:
                result = int(reads_value_list[second_max_index])/int(reads_value_list[max_index])
            else:
                result = int(reads_value_list[max_index])/int(reads_value_list[second_max_index])
            reads_value_ratio.append(result)
            reads_value_ratio_dic[i[0]]=str(result)+'\t'+reads
            reads_value_list=[]
    return reads_value_list, reads_value_ratio

def get_count_reads_ratio(input_file):
    global value_list
    global max_index
    global second_max_index
    global result
    global max_value
    global second_max_value
    value_list = []
    value_ratio1 = []
    value_ratio1_dic = {}
    value_ratio2 = []
    value_ratio2_dic = {}
    value_ratio3 = []
    value_ratio3_dic = {}
    value_ratio3_inverse = []
    value_ratio3_inverse_dic = {}
    value_ratio4 = []
    value_ratio4_dic = {}
    f= open(input_file,'r')
    for i in f:
        j=i
        #print(j)
        i=i.split('\t')
        reads=i[1]+'\t'+i[2]+'\t'+i[3]+'\t'+i[4].strip('\n')
        #print(reads)
        value_list.extend((i[1],i[2],i[3],i[4].strip('\n')))
        value_list = list(map(lambda x: int(x), value_list))
        #print(value_list)
        max_index = value_list.index(max(value_list))
        second_max_index = value_list.index(sorted(value_list)[-2])
        max_value = value_list[max_index]
        second_max_value = value_list[second_max_index]
        #print(str(max_value)+'\t'+str(second_max_value))
        if second_max_value != 0:
            if max_index > second_max_index:
                result = int(value_list[second_max_index])/int(value_list[max_index])
            else:
                result = int(value_list[max_index])/int(value_list[second_max_index])
            value_ratio1.append(result)
            value_ratio1_dic[i[0]]=str(result)+'\t'+reads
            if value_list[second_max_index] <= 2:
                #print(result)
                value_ratio3_inverse.append(result)
                value_ratio3_inverse_dic[i[0]]=str(result)+'\t'+reads
            else:
                value_ratio3.append(result)
                value_ratio3_dic[i[0]]=str(result)+'\t'+reads
            #print(i[0]+'\t'+str(result))
            value_list=[]
    #print(len(value_ratio1))
    return value_ratio1,value_ratio1_dic,value_ratio3,value_ratio3_inverse,value_ratio3_dic,value_ratio3_inverse_dic,value_ratio2,value_ratio2_dic,value_ratio4,value_ratio4_dic


def get_final_four_thresholds_snps(reads_value_ratio ,a_value_ratio1, a_value_ratio1_dic, a_value_ratio2, a_value_ratio2_dic, a_value_ratio3, a_value_ratio3_dic, a_value_ratio3_inverse, a_value_ratio3_inverse_dic, a_value_ratio4, a_value_ratio4_dic):
    median = np.percentile(a_value_ratio1, 50)
    lower = np.percentile(a_value_ratio1, 10)
    upper = np.percentile(a_value_ratio1, 90)
    for a in a_value_ratio1_dic.keys():
        if float(a_value_ratio1_dic[a].split('\t')[0]) > lower and float(a_value_ratio1_dic[a].split('\t')[0]) < upper :
            a_value_ratio2.append(float(a_value_ratio1_dic[a].split('\t')[0]))
            a_value_ratio2_dic[a]=a_value_ratio1_dic[a]
        else:
            pass
    for a in a_value_ratio3_dic.keys():
        if float(a_value_ratio3_dic[a].split('\t')[0]) > lower and float(a_value_ratio3_dic[a].split('\t')[0]) < upper :
            a_value_ratio4.append(float(a_value_ratio3_dic[a].split('\t')[0]))
            a_value_ratio4_dic[a]=a_value_ratio3_dic[a]
        else:
            pass
    
    # dfo=open(tmp_dir+'/each_snp_coverage_barcodes_filtered.txt','w')
    dfo_inter1 = open(tmp_dir + '/each_snp_coverage_barcodes_ratio.txt', 'w')
    dfo_inter2 = open(tmp_dir + '/each_snp_coverage_barcodes_filtered_by_ratio_ratio.txt', 'w')
    dfo_inter3 = open(tmp_dir + '/each_snp_coverage_barcodes_filtered_by_reads_ratio.txt', 'w')
    dfo_inter3_inverse = open(tmp_dir + '/each_snp_coverage_barcodes_unmeet_by_reads_ratio.txt', 'w')
    dfo_inter4 = open(tmp_dir + '/each_snp_coverage_barcodes_filtered_by_read_and_ratio_ratio.txt', 'w')

    for i in a_value_ratio1_dic.keys():
        #print(str(value_ratio1_dic[i]))
        dfo_inter1.write(i+'\t'+str(a_value_ratio1_dic[i])+'\n')
    for i in a_value_ratio2_dic.keys():
        #print(str(value_ratio2_dic[i]))
        dfo_inter2.write(i+'\t'+str(a_value_ratio2_dic[i])+'\n')
    for i in a_value_ratio3_dic.keys():
        #print(str(value_ratio3_dic[i]))
        dfo_inter3.write(i+'\t'+str(a_value_ratio3_dic[i])+'\n')
    for i in a_value_ratio4_dic.keys():
        #print(str(value_ratio4_dic[i]))
        dfo_inter4.write(i+'\t'+str(a_value_ratio4_dic[i])+'\n')

    ### visulization
    log2_unfiltered_list = [math.log2(x) for x in a_value_ratio1]
    log2_filtered_ratio2_list = [math.log2(x) for x in a_value_ratio2]
    log2_filtered_ratio3_list = [math.log2(x) for x in a_value_ratio3]
    log2_filtered_ratio4_list = [math.log2(x) for x in a_value_ratio4]
    log2_unfiltered_reads_list = [math.log2(x) for x in reads_value_ratio]

    ###### unfold
    ### unfiltered reads
    #n, bins, patches = plt.hist(value_ratio3, bins=40)
    plt.hist(a_value_ratio3,density=True,bins=40)
    #sns.kdeplot(value_ratio3,bw_method=0.05, label='Manual')
    plt.xlabel('ratio Value[reads]')
    plt.ylabel('Density')
    plt.title('Distribution of unfiltered snp\'s coverage by reads number')
    plt.savefig(tmp_dir+'/unfiltered_reads_data_distribution.pdf')
    plt.clf()

    ## unfiltered barcodes
    #n, bins, patches = plt.hist(value_ratio1, bins=40)
    plt.hist(a_value_ratio1,density=True,bins=40)
    #sns.kdeplot(value_ratio1,bw_method=0.05, label='Manual')
    plt.xlabel('ratio Value')
    plt.ylabel('Density')
    plt.title('Distribution of unfiltered snp\'s coverage')
    plt.savefig(tmp_dir+'/unfiltered_data_distribution.pdf')
    plt.clf()

    ## filtered by ratio
    plt.hist(a_value_ratio2,density=True,bins=40)
    #sns.kdeplot(value_ratio2,bw_method=0.05, label='Manual')
    plt.xlabel('ratio Value')
    plt.ylabel('Density')
    plt.title('Distribution of filtered snp\'s coverage by ratio')
    plt.savefig(tmp_dir+'/filtered_data_by_ratio_distribution.pdf')
    plt.clf()

    ## filtered by ratio and reads
    plt.hist(a_value_ratio4,density=True,bins=40)
    #sns.kdeplot(value_ratio4,bw_method=0.05, label='Manual')
    plt.xlabel('ratio Value')
    plt.ylabel('Density')
    plt.title('Distribution of filtered snp\'s coverage by ratio and reads')
    plt.savefig(tmp_dir+'/filtered_data_by_ratio_and_reads_distribution.pdf')
    plt.clf()

    ######
    ### log2
    ### unfiltered reads
    #n, bins, patches = plt.hist(log2_unfiltered_reads_list, bins=40)
    plt.hist(log2_unfiltered_reads_list,density=True,bins=40)
    #sns.kdeplot(log2_unfiltered_reads_list,bw_method=0.05, label='Manual')
    plt.xlabel('ratio Value[reads]')
    plt.ylabel('Density')
    plt.title('Distribution of unfiltered snp\'s coverage by reads number')
    plt.savefig(tmp_dir+'/unfiltered_reads_data_log2_distribution.pdf')
    plt.clf()

    ## unfiltered barcodes
    #n, bins, patches = plt.hist(log2_unfiltered_list, bins=40)
    plt.hist(log2_unfiltered_list,density=True,bins=40)
    #sns.kdeplot(log2_unfiltered_list,bw_method=0.05, label='Manual')
    plt.xlabel('ratio Value')
    plt.ylabel('Density')
    plt.title('Distribution of unfiltered snp\'s coverage')
    plt.savefig(tmp_dir+'/unfiltered_data_log2_distribution.pdf')
    plt.clf()

    ## filtered by ratio
    plt.hist(log2_filtered_ratio2_list,density=True,bins=40)
    #sns.kdeplot(log2_filtered_ratio2_list,bw_method=0.05, label='Manual')
    plt.xlabel('ratio Value')
    plt.ylabel('Density')
    plt.title('Distribution of filtered snp\'s coverage by ratio')
    plt.savefig(tmp_dir+'/filtered_data_by_ratio_log2_distribution.pdf')
    plt.clf()

    ## filtered by ratio and reads
    plt.hist(log2_filtered_ratio4_list,density=True,bins=40)
    #sns.kdeplot(log2_filtered_ratio4_list,bw_method=0.05, label='Manual')
    plt.xlabel('ratio Value')
    plt.ylabel('Density')
    plt.title('Distribution of filtered snp\'s coverage by ratio and reads')
    plt.savefig(tmp_dir+'/filtered_data_by_ratio_and_reads_log2_distribution.pdf')
    plt.clf()

    print('number of unfiltered snp is ',str(len(a_value_ratio1)))
    print('get distribution of unfiltered snp done!')
    print('number of filtered snp by ratio is ',str(len(a_value_ratio2)))
    print('get distribution of filtered snp by ratio done!')
    print('number of filtered snp by reads is ',str(len(a_value_ratio3)))
    print('get distribution of filtered snp by reads done!')
    print('number of filtered snp by ratio and reads is ',str(len(a_value_ratio4)))
    print('get distribution of filtered snp by ratio and reads done!')


################################################# step6: get result_matrix1.csv (merged snp-barcode matrix)
# define function that can get merged all snp base info of each barcode (snp-barcode matrix1)
# sgtxt_list = []
# all_path(inter_base_dir,'merged_base_count.txt')
# barcodes_ACGT_count_list = sgtxt_list
def get_snp_barcode_matrix1_dir(z):
    global inter_base_dir
    global sgtxt_list
    final_path = []
    count_files = all_path(inter_base_dir,'merged_base_count.txt')
    for count_file in count_files:
        tmp_final_path = os.path.join(inter_base_dir, count_file)
        final_path.append(tmp_final_path)
    #print(len(final_path))
    return final_path

                             # final path
def get_snp_barcode_matrix1(path):
    df_list=[]
    dict_res = {}
    barcode_list=[]
    snp_list=[]
    snp_list1=[]
    n=0
    #print(len(path))
    for i,file in enumerate(path):
        #print(i)
        n+=1
        snp_list=[]
        name=os.path.basename(file)
        barcode=name[0:18]
        barcode_list.append(barcode)
        df = pd.read_csv(file,sep='\t', header=None, names=['Chromosome:Position', 'A;C;G;T'])
        df_list.append(df)
        snp_list=df['Chromosome:Position'].values
        dict_res[name]=df[['A;C;G;T']]
        if i == 1:
            snp_list1=snp_list
    #print(dict_res.values())
    #return df_list, snp_list1, barcode_list
    column = 'A;C;G;T'
    #dfs = list(df_list)
    #df_list = [df[~df.index.duplicated()] for df in df_list]
    dfs = df_list
    dfs_flat = [df.set_index(df.index.to_flat_index()) for df in dfs]
    result = pd.concat([df[column] for df in dfs_flat], axis=1,ignore_index=False)
    data = result

    data.columns=barcode_list
    data.index=snp_list1
    #print('create barcode-snp matrix1 done!')

    data_t = data.T
    #data_t2 = data_t1.replace('NaN', '.;.;.;.')
    data_t2 = data_t.fillna(".;.;.;.")
    #print(data_t2)
    data_t2.to_csv(result_dir+'/result_matrix1.csv',sep='\t',mode='w',header=True,index=True)
    print(str('din for matrix1 is ')+str(data_t2.shape))
    print('get filtered barcode-snp matrix1 done!')
    
    print(str('length for df_list is ')+str(len(df_list)))
    return data_t2


################################################### step7: get result_matrix2.csv, result_matrix3.csv, and result_matrix4.csv
####  create matrix1_supple after base filtering (matrix2, matrix3, and matrix4)
# matrix1:filtering by only sequence charateristics (upstream and downstream 10bp)
# matrix2:filtering by only reads(second max base counts >2)
# matrix3:filtering by only ratio(0.1-0.9)
# matrix4:filtering by both reads(second max base base counts >2) and ratio(0.1-0.9)
# 如上上个矩阵(matrix2,matrix3,matrix4)分别对应输出的filter_df1,filter_df2,filter_df3
def get_snp_barcode_matrix2_matrix3_matrix4(matrix1):
    # chrX;20282118   2.6666666666666665      0       0       32      12
    df_filter2 = pd.read_csv(tmp_dir+'/each_snp_coverage_barcodes_filtered_by_reads_ratio.txt',sep='\t',header=None, names=['Chromosome;Position', 'ratio','A','C','G','T'])
    df_filter3 = pd.read_csv(tmp_dir+'/each_snp_coverage_barcodes_filtered_by_ratio_ratio.txt',sep='\t',header=None, names=['Chromosome;Position', 'ratio','A','C','G','T'])
    df_filter4 = pd.read_csv(tmp_dir+'/each_snp_coverage_barcodes_filtered_by_read_and_ratio_ratio.txt',sep='\t', header=None, names=['Chromosome;Position', 'ratio','A','C','G','T'])
    content_list2 = df_filter2['Chromosome;Position'].tolist()
    content_list3 = df_filter3['Chromosome;Position'].tolist()
    content_list4 = df_filter4['Chromosome;Position'].tolist()
    #print(content_list4)

    df = matrix1
    df_list = df.columns.tolist()

    filter_list2 = list(set(df_list) & set(content_list2))
    filter_list3 = list(set(df_list) & set(content_list3))
    filter_list4 = list(set(df_list) & set(content_list4))
    #print(len(filter_list2))
    #print(len(filter_list3))
    #print(len(filter_list4))

    # matrix2
    filtered_df2 = df[filter_list2]
    # matrix3
    filtered_df3 = df[filter_list3]
    # matrix4
    filtered_df4 = df[filter_list4]

    row_names = df.index.tolist()
    filtered_df2.index = row_names
    filtered_df3.index = row_names
    filtered_df4.index = row_names
    print('dims for barcode-snp matrix2 is ',str(filtered_df2.shape))
    print('dims for barcode-snp matrix3 is ',str(filtered_df3.shape))
    print('dims for barcode-snp matrix4 is ',str(filtered_df4.shape))

    filtered_df2.to_csv(result_dir+'/result_matrix2.csv',sep='\t',mode='w',header=True,index=True)
    filtered_df3.to_csv(result_dir+'/result_matrix3.csv',sep='\t',mode='w',header=True,index=True)
    filtered_df4.to_csv(result_dir+'/result_matrix4.csv',sep='\t',mode='w',header=True,index=True)
    print('get filtered barcode-snp matrix2, matrix3, matrix4 done!')


############################################# step8 : delete temporary files 
def delete_files_in_directory(directory, file_endfix):
    for filename in os.listdir(directory):
        if filename.endswith(file_endfix):
            file_path = os.path.join(directory, filename)
            os.remove(file_path)
            print(f"Deleted file: {file_path}")




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
    all_path(inter_base_dir, 'bam')
    for i in sgtxt_list:
        small_bam_file = inter_base_dir + '/' + i
        subprocess.run([samtools, "index", small_bam_file])
    print('create index for all the barcodes, done!')


    ############################# step2: get ACGT count of all analysed SNPs for each barcode
    sgtxt_list = []
    all_path(inter_base_dir, '.bam')
    bam_files = sgtxt_list
    genome = args.genome
    bed_file = tmp_dir + '/filtered_freebayes_snps_manual_info.txt'
    bed_file2=tmp_dir+'/filtered_freebayes_snps.bed'
    get_each_snp_ACGT_count_for_each_barcode(bam_files, bed_file, genome)


    ############################## step3: get merged ACGT count of all analysed SNPs for each barcode
    manual_file=tmp_dir+'/filtered_freebayes_snps_manual_info.txt'
    data = get_merged_snps_ACGT_count_file_for_each_barcode(manual_file)


    ############################## step4: get snp_information
    snp_loci_list = []
    snp_loci_dic = {}
    dic1 = {}
    dic2 = {}
    snp_loci_list1 = []
    snp_loci_list2 = []

    input_file = tmp_dir + '/filtered_freebayes_snps_manual_info.txt'
    snp_loci_list, snp_loci_dic = define_filtered_snp_loci(input_file)
    dic1 = snp_loci_dic
    snp_loci_list1 = snp_loci_list
    print(str('length of dic1 is :')+ str(len(dic1))+str('!'))

    snp_loci_list = []
    snp_loci_dic = {}
    input_file = tmp_dir + '/filtered_freebayes_snps.bed'
    snp_loci_list, snp_loci_dic = define_filtered_snp_loci(input_file)
    dic2 = snp_loci_dic
    snp_loci_list2 = snp_loci_list
    print(str('length of dic2 is :')+ str(len(dic2))+str('!'))

    # # filtered by snp and sequence by manual methods,reads and ratio
    common_dir = inter_base_dir
    sgtxt_list = all_path(inter_base_dir, '-1_base_count.txt')
    print(len(sgtxt_list))

    each_snp_info_dic,each_snp_info_dic1,snp_loci_4base_barcodes_list_dic = get_single_snp_coverage(sgtxt_list, dic1, snp_loci_list1)
    
    foa = open(tmp_dir + '/each_snp_coverage_barcodes.txt', 'w')
    for key in each_snp_info_dic.keys():
        string = '\t'.join('%s' % id for id in each_snp_info_dic[key])
        foa.write(key + '\t' + string + '\n')
    print('count each filtered snp coverages in filtered barcodes done!')

    foa1 = open(tmp_dir + '/each_snp_coverage_reads.txt', 'w')
    for key in each_snp_info_dic1.keys():
        string1 = '\t'.join('%s' % id for id in each_snp_info_dic1[key])
        foa1.write(key + '\t' + string1 + '\n')
    print('count each filtered snp reads coverages in filtered barcodes done!')


    ############################## step5: get two main bases's distribution of each snp (unfiltered)
    input_file = tmp_dir + '/filtered_freebayes_snps_manual_info.txt'
    snp_sequence_dic = get_final_snp_dic(input_file)

    ### up and down20bp sequence of snp after filtering snp by manual parameters
    # dfo_inter5=open(tmp_dir+'/filtered_freebayes_snps_manual_info_by_ratio_and_reads.txt','w')
    # chrX;299510     31      0       3       0
    input_file2 = tmp_dir + '/each_snp_coverage_reads.txt'
    reads_value_list, reads_value_ratio = get_unfiltered_reads(input_file2)

    input_file = tmp_dir + '/each_snp_coverage_barcodes.txt'
    #value_ratio1,value_ratio1_dic,value_ratio2,value_ratio2_dic,value_ratio3,value_ratio3_inverse,value_ratio3_dic,value_ratio3_inverse_dic,value_ratio4,value_ratio4_dic = get_count_reads_ratio(input_file)
    value_ratio1,value_ratio1_dic,value_ratio3,value_ratio3_inverse,value_ratio3_dic,value_ratio3_inverse_dic,value_ratio2,value_ratio2_dic,value_ratio4,value_ratio4_dic = get_count_reads_ratio(input_file)
    #print(len(value_ratio1))

    get_final_four_thresholds_snps(reads_value_ratio, value_ratio1, value_ratio1_dic, value_ratio2, value_ratio2_dic, value_ratio3, value_ratio3_dic, value_ratio3_inverse, value_ratio3_inverse_dic, value_ratio4, value_ratio4_dic)


    ############################## step6: get result_matrix1.csv (merged snp-barcode matrix)
    final_path = get_snp_barcode_matrix1_dir(1)

    matrix1 = get_snp_barcode_matrix1(final_path)


    ############################## step7: get result_matrix2.csv, result_matrix3.csv, and result_matrix4.csv
    get_snp_barcode_matrix2_matrix3_matrix4(matrix1)
    

    ############################## step8: delete temporary files
    directory_path = inter_base_dir
    file_endfix_to_delete = "bam"
    delete_files_in_directory(directory_path, file_endfix_to_delete)
    file_endfix_to_delete = "bam.bai"
    delete_files_in_directory(directory_path, file_endfix_to_delete)
    file_endfix_to_delete = "base_count.txt"
    delete_files_in_directory(directory_path, file_endfix_to_delete)
    shutil.rmtree(inter_base_dir)

    directory_path = tmp_dir
    file_endfix_to_delete = "txt"
    delete_files_in_directory(directory_path, file_endfix_to_delete)
    file_endfix_to_delete = "vcf"
    delete_files_in_directory(directory_path, file_endfix_to_delete)
    file_endfix_to_delete = "bed"
    delete_files_in_directory(directory_path, file_endfix_to_delete)
    file_endfix_to_delete = "bam"
    delete_files_in_directory(directory_path, file_endfix_to_delete)
    file_endfix_to_delete = "bam.bai"
    delete_files_in_directory(directory_path, file_endfix_to_delete)

    print('preprocess data,done!')

x=1
runpipeline(x)



