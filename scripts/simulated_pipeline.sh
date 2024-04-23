### 1：Download the source data GSE203275 and run CellRanger. 

### 2: Add the MD labels to align results.
nohup samtools calmd -@ 36 -bAr possorted_genome_bam.bam /path_to_fa/genome.fa > possorted_baq_ne_genome_bam.bam &

### 3：Change the reads to reference, and then change to the SNP of the HG000405 in the 1000 genome trio to distinguish the origin of the parent.
nohup python change_to_ref_to_405_v2,py &
