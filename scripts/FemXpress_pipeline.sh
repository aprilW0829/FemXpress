#!usr/bin/bash

# step1:
# data process
python FemXpress_1.py possorted_genome_bam.bam mm10.fa metadata.txt mm10.rmsk.txt 

# step2:
# FemXpress inference clustering and escape genes
python FemXpress_2.py result_matrix4.txt > nohup_real.txt
