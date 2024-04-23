#!usr/bin/bash

# step1:
# data preprocess
python FemXpress_1.py -b possorted_genome_bam.bam -g genome.fa -e meta.txt -r rmsk.txt > nohup_real1.txt

# step2:
# FemXpress inference clustering and escape genes
python FemXpress_2.py -m result_matrix4.csv -a /data2/wangxin/database/genome/gencode.vM25.chr_patch_hapl_scaff.basic.annotation.gtf > nohup_real2.txt
