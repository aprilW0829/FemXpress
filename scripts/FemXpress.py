#!usr/bin/python3
# -*- coding: utf-8 -*-
""""
usage: FemXpress.py [-h] [-b BAM] [-g GENOME] [-e META] [-r RMSK] [-n SAMPLE] [-a ANNOT] [-p PARENTS_GENOTYPES] [-c MINIMUM_COUNTS] [-v]

FemXpress anlysis, Example:python FemXpress.py -b possorted_genome_bam.bam -g genome.fa -e meta.txt -r rmsk.txt -n sample_name -a gene_annotation.gtf -p
parents_genotypes.csv -c 2

optional arguments:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     input bam file from CellRanger output
  -g GENOME, --genome GENOME
                        input genome
  -e META, --meta META  input meta information of barcodes
  -r RMSK, --rmsk RMSK  rmsk file
  -n SAMPLE, --sample SAMPLE
                        sample name
  -a ANNOT, --annot ANNOT
                        gene annotation file
  -p PARENTS_GENOTYPES, --parents_genotypes PARENTS_GENOTYPES
                        parents genotypes file
  -c MINIMUM_COUNTS, --minimum_counts MINIMUM_COUNTS
                        minimum counts supporting second largest base, default: 2
  -v, --verbose         print verbose output
"""

import os
import sys
import threading
from time import time
import argparse

def run_script(sample=None, bam=None, genome=None, meta=None, rmsk=None, annot=None, parents_genotypes=None, minimum_counts=None, SNP_qc=None):

    start_time = time()
    
    basename = sample
    os.system("mkdir FemXpress/log")
    '''
    print("######## 01 FemXpress analysis, data preprocess, step1... ########")
    cmd = "python FemXpress_preprocess1.py -b {} -g {} -e {} -r {} -n {} -c {}".format(bam, genome, meta, rmsk, sample, minimum_counts)
    print(cmd)
    os.system(cmd)
    print("FemXpress preprocess step1 Done in {}s".format(round(time()-start_time, 3)))

    print("######## 02 FemXpress analysis, data preprocess, step2... ########")
    cmd = "python FemXpress_preprocess2.py -b {} -g {} -e {} -r {} -n {} -c {}".format(bam, genome, meta, rmsk, sample, minimum_counts)
    print(cmd)
    os.system(cmd)
    print("FemXpress preprocess step2 Done in {}s".format(round(time()-start_time, 3)))
    '''
    if SNP_qc:
        matrix = './FemXpress/'+sample+'/result/result_matrix4.csv'
    else:
        matrix = './FemXpress/'+sample+'result/result_matrix1.csv'
    print("######## 02 FemXpress analysis... ########")
    cmd = "python FemXpress_process.py -m {} -a {} -n {} -p {}".format(matrix, annot, sample, parents_genotypes)
    print(cmd)
    os.system(cmd)

    print("All done, please check log files!")
    print("FemXpress Done in {}s".format(round(time()-start_time, 3)))


def main():
    parser = argparse.ArgumentParser(description='FemXpress anlysis,\nExample:python FemXpress.py -b possorted_genome_bam.bam -g genome.fa -e meta.txt -r rmsk.txt -n sample_name -a gene_annotation.gtf -p parents_genotypes.csv -c 2')
    parser.add_argument("-b","--bam", type=str,help="input bam file from CellRanger output")
    parser.add_argument("-g","--genome",type=str,help="input genome")
    parser.add_argument("-e","--meta",type=str,help="input meta information of barcodes")
    parser.add_argument("-r","--rmsk",type=str,help="rmsk file")
    parser.add_argument("-n","--sample",type=str,help="sample name")
    parser.add_argument("-a","--annot",type=str,help="gene annotation file")
    parser.add_argument("-p","--parents_genotypes",type=str,help="parents genotypes file")
    parser.add_argument("-c","--minimum_counts", type=int,help="minimum counts supporting second largest base, default: 2")
    parser.add_argument("-q","--SNP_qc", action="store_true", dest="SNP_qc", default=True, help="remove SNPs with 10% highest and lowest ratios")
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=False, help="print verbose output")
    args=parser.parse_args()

    print('[FemXpress] genome: ' + str(args.genome))
    print('[FemXpress] sample: ' + str(args.sample))
    print('[FemXpress] bam: ' + str(args.bam))
    print('[FemXpress] meta: ' + str(args.meta))
    print('[FemXpress] rmsk: ' + str(args.rmsk))
    print('[FemXpress] annot: ' + str(args.annot))
    print('[FemXpress] parents_genotypes: ' + str(args.parents_genotypes))
    print('[FemXpress] minimum_counts: ' + str(args.minimum_counts))

    try:
        run_script(args.sample, args.bam, args.genome, args.meta, args.rmsk, args.annot, args.parents_genotypes, args.minimum_counts)
    except KeyboardInterrupt:
        sys.stderr.write("See you~ :)\n")
        sys.exit(0)

if __name__ == "__main__":
    main()
