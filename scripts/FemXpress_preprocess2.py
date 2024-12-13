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
print("[FenXpress preprocess] file number's limit:", soft_limit, hard_limit)
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


############################################ step3: get merged ACGT count of all analysed SNPs for each barcode
# define dunction that can get equal snps for each barcode (using whole snps, can be ,;,;,;,)
def get_merged_snps_ACGT_count_file_for_each_barcode(manual_file,count_files):
    
    data = [] 
    positions = []
    global inter_base_dir
    global tmp_dir

    #with open(tmp_dir+'/filtered_freebayes_snps_manual_info.txt', 'r') as bed_file:
    with open(manual_file,'r') as bed_file:
        for line in bed_file:
            chrom, start, end, up_down_20bp_sequence = line.strip('\n').split('\t')
            a=str(chrom)+';'+str(start)
            data.append(a)
            positions.append((chrom,int(start)))
    data=list(set(data))
    #print(data[:6])
    print('[FemXpress preprocess] number of snp analysed is '+str(len(data)))
    #return data

    columns_names=['loci','A;C;G;T','ref','total_count']
    for count_file in count_files:
        df = pd.read_csv(inter_base_dir+'/'+count_file,names=columns_names,sep='\t')
        fw = open(inter_base_dir+'/'+count_file[:18]+'_merged_base_count.txt','w')
        df = df.drop_duplicates(subset=['loci'])
        search_col = 'loci'
        output_col = 'A;C;G;T'
        for search_value in data:
            #print(search_value)
            if search_value in df[search_col].values:
                row=df.loc[df[search_col]==search_value]
                #result1 = np.array2string(row[['loci','A;C;G;T']].values,separator='\t', formatter={'int_kind':lambda x: "%.2f" % x})
                result1 = np.array2string(row[['loci','A;C;G;T']].values,separator='\t')
                result1 = result1.replace('[','').replace(']','').replace("'","").strip("'")
                #print(result1)
                #print('there is result1!')
                #if result1.split('\t')[1]=='0;0;0;0':
                if result1.split('\t')[1] == '0;0;0;0':
                    fw.write(result1.split('\t')[0]+ '\t.;.;.;.\n')
                elif result1.split('\t')[1]==' ':
                    fw.write(result1.split('\t')[0]+ '\t.;.;.;.\n')
                else:
                    #print(result1)
                    #print('there is result1!')
                    fw.write(result1+ '\n')
            else:
                fw.write(search_value+'\t.;.;.;.\n')
       # print(count_file, 'done!')

    print('get merged_base_count for each barcode done!')
    return data


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
            if value_list[second_max_index] <= args.minimum_counts:
                #print(result)
                value_ratio3_inverse.append(result)
                value_ratio3_inverse_dic[i[0]]=str(result)+'\t'+reads
            else:
                value_ratio3.append(result)
                value_ratio3_dic[i[0]]=str(result)+'\t'+reads
            #print(i[0]+'\t'+str(result))
            value_list=[]
    #print(len(value_ratio1))
    #print(str(len(value_ratio1)), str(len(value_ratio1_dic)), str(len(value_ratio3)), str(len(value_ratio3_inverse)), str(len(value_ratio3_dic)), str(len(value_ratio2)), str(len(value_ratio2_dic)), str(len(value_ratio4)), str(len(value_ratio4_dic)))
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

    print('[FemXpress preprocess] number of unfiltered snp is ',str(len(a_value_ratio1)))
    print('[FemXpress preprocess] get distribution of unfiltered snp done!')
    print('[FemXpress preprocess] number of filtered snp by ratio is ',str(len(a_value_ratio2)))
    print('[FemXpress preprocess] get distribution of filtered snp by ratio done!')
    print('[FemXpress preprocess] number of filtered snp by reads is ',str(len(a_value_ratio3)))
    print('[FemXpress preprocess] get distribution of filtered snp by reads done!')
    print('[FemXpress preprocess] number of filtered snp by ratio and reads is ',str(len(a_value_ratio4)))
    print('[FemXpress preprocess] get distribution of filtered snp by ratio and reads done!')


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
    for column in data_t2.columns:
        if (data_t2[column] == '.;.;.;.').all():
            data_t2 = data_t2.drop(column, axis=1)

    data_t2.to_csv(result_dir+'/result_matrix1.csv',sep='\t',mode='w',header=True,index=True)
    print('[FemXpress preprocess] '+str('shape for matrix1 is ')+str(data_t2.shape))
    print('[FemXpress preprocess] get filtered barcode-snp matrix1 done!')
    
    #print(str('length for df_list is ')+str(len(df_list)))
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
    #df = df[df['column_name'] != '.;.;.;']

    print('[FemXpress preprocess] shape for barcode-snp matrix2 is ',str(filtered_df2.shape))
    print('[FemXpress preprocess] shape for barcode-snp matrix3 is ',str(filtered_df3.shape))
    print('[FemXpress preprocess] shape for barcode-snp matrix4 is ',str(filtered_df4.shape))

    filtered_df2.to_csv(result_dir+'/result_matrix2.csv',sep='\t',mode='w',header=True,index=True)
    filtered_df3.to_csv(result_dir+'/result_matrix3.csv',sep='\t',mode='w',header=True,index=True)
    filtered_df4.to_csv(result_dir+'/result_matrix4.csv',sep='\t',mode='w',header=True,index=True)
    print('[FemXpress preprocess] get filtered barcode-snp matrix2, matrix3, matrix4 done!')


############################################# step8 : delete temporary files 
def delete_files_in_directory(directory, file_endfix):
    for filename in os.listdir(directory):
        if filename.endswith(file_endfix):
            file_path = os.path.join(directory, filename)
            os.remove(file_path)
            #print(f"Deleted file: {file_path}")




# define function that run pipeline
def runpipeline(a):
    global tmp_dir
    global inter_base_dir
    global result_dir
    global samtools
    global trimmed_barcode_list
    global sgtxt_list 

    ############################## step3: get merged ACGT count of all analysed SNPs for each barcode
    manual_file=tmp_dir+'/filtered_freebayes_snps_manual_info.txt'
    sgtxt_list=[]
    sgtxt_list = all_path(inter_base_dir,'-1_base_count.txt')
    count_files=sgtxt_list
    print('[FemXpress preprocess] length of analysed base_count files is '+str(len(count_files)))
    data = get_merged_snps_ACGT_count_file_for_each_barcode(manual_file,count_files)


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
    print('[FemXpress preprocess] ' +str('length of dic1 is :')+ str(len(dic1))+str('!'))

    snp_loci_list = []
    snp_loci_dic = {}
    input_file = tmp_dir + '/filtered_freebayes_snps.bed'
    snp_loci_list, snp_loci_dic = define_filtered_snp_loci(input_file)
    dic2 = snp_loci_dic
    snp_loci_list2 = snp_loci_list
    print('[FemXpress preprocess] ' +str('length of dic2 is :')+ str(len(dic2))+str('!'))

    # # filtered by snp and sequence by manual methods,reads and ratio
    common_dir = inter_base_dir
    sgtxt_list = all_path(inter_base_dir, '-1_base_count.txt')
    #print(len(sgtxt_list))

    each_snp_info_dic,each_snp_info_dic1,snp_loci_4base_barcodes_list_dic = get_single_snp_coverage(sgtxt_list, dic1, snp_loci_list1)
    
    foa = open(tmp_dir + '/each_snp_coverage_barcodes.txt', 'w')
    for key in each_snp_info_dic.keys():
        string = '\t'.join('%s' % id for id in each_snp_info_dic[key])
        foa.write(key + '\t' + string + '\n')
    print('[FemXpress preprocess] count each filtered snp coverages in filtered barcodes done!')

    foa1 = open(tmp_dir + '/each_snp_coverage_reads.txt', 'w')
    for key in each_snp_info_dic1.keys():
        string1 = '\t'.join('%s' % id for id in each_snp_info_dic1[key])
        foa1.write(key + '\t' + string1 + '\n')
    print('[FemXpress preprocess] count each filtered snp reads coverages in filtered barcodes done!')


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

    print('[FemXpress preprocess] preprocess data, step2 ,done!')

x=1
runpipeline(x)
