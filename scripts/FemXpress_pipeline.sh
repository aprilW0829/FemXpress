#!usr/bin/bash

# step1:
# data preprocess
python FemXpress_1.py -b possorted_genome_bam.bam -g genome.fa -e meta.txt -r rmsk.txt

# step2:
# FemXpress inference clustering and escape genes
python FemXpress_2.py result_matrix4.txt > nohup_real.txt
