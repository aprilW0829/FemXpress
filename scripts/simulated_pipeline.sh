### 1：下载源数据GSE203275，运行cellranger 

### 2: 把比对结果加上MD标签 
# nohup samtools calmd -@ 36 -bAr possorted_genome_bam.bam /path_to_fa/genome.fa > possorted_baq_ne_genome_bam.bam &

### 3：把reads修改为reference,再修改为1000 genome trio中的HG000405的能分辨出父母来源的SNP
nohup python change_to_ref_to_405_v2,py &
