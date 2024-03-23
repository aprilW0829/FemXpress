# FemXpress
FemXpress--tool that classifies the source of X chromosome inactivation of each female cell 

## User manual and guide
### Overview
FemXpress is a powerful tool designed for a comprehensive analysis of XCI status in female scRNA-Seq datasets. It has the following main functions: 1. to infer the source of X chromosome inactivation (paternal or maternal) for each cell, 2. to infer new potential escape genes, and 3. to infer the X-to-auto expression dose in different cell types.

To avoid relying on parental SNP information, we performed de novo discovery of SNP that exhibit parental origin differences using the alignment results of the pseudo-bulk scRNA-Seq sample. In order to minimize the impact of sequencing and analysis errors on SNP calling, we applied strict filtering criteria, including: 1) a minimum coverage of at least 2 reads supporting either of the two genotypes of the SNP, 2) a minimum number of 2 cells with each genotype, 3) exclusion of loci originating from simple repeat regions, and 4) sorting the remaining potential SNPs based on the ratio of reads for each genotype at the locus and removing the top and bottom 10% of SNPs with the highest and lowest ratios, respectively, as they are likely introduced by sequencing errors in high coverage regions (Fig. 1a-c).
Subsequently, linkage analysis was performed based on the identified high-confidence SNPs. We counted the genotypes of two SNPs that appear simultaneously in each cell, considering them linked genotypes only if both genotypes are present in at least 80% of cells (Fig. 1d). Next, we conducted an extended linkage scoring by combining linked genotypes into longer chains of linked loci. The pair of haplotypes with the highest number of linked loci was selected as the genotypes of the parental source X chromosome (Fig. 1e). Then, when the haplotypes of the parents were phased, the origin of the X chromosome for each cell was determined by a voting method. Each haplotype received one vote if the cell carried an allele base from that haplotype, and the cell was classified to have that haplotype if the votes were higher than those of the other haplotype. Finally, in the cells that were already classified with the same parental source of the inactivated X chromosome, we examined the frequency of concurrent parental origin genotypes at each SNP locus to identify potential XCI escaping gene derived SNPs(Fig. 1f).
![Image Alt Text](image_file_path)


## Citing FemXpress
A preprint FemXpress: Systematic Analysis of X Chromosome Inactivation Heterogeneity in Female Single-Cell RNA-Seq Samples provides an overview of the statistical models used in FemXpress. Ww ask that you cite this paper if you use FemXpress in work that leads to publication.This preprint is used for documentation and citation.

## Install
To install FemXpress,make sure you need install some dependances,if you need more details on the dependences,look at the FemXpress-env.yaml file.
set up conda environment for FemXpress
conda env create -f FemXpress-env.yaml

## install FemXpress from zip package
conda activate FemXpress-env
wget http://www.example.com/FemXpress.zip
gunzip FemXpress.zip

## Usage
The current version of FemXpress was only developed based on single-cell sequencing data from 10×Genomics. It consists of two running steps:
Step 1 : to obtain pre-processed barcode-SNP matrix based on the alignment bam file from cellranger’s output.
According to how it works, it needs to have the following input files: BAM file for CellRanger alignment output, cell tag file (which can be meta information file that has been processed by Seurat)(possorted_genome_bam.bam), reference genome sequence file(.fa), and rmsk file.For instance:

Python FemXpress_1.py possorted_genome_bam.bam genome.fa meta.txt rmsk.txt

...will produce four barcode-SNP matrix files in the subdirectory ./FemXpress/result, result_matrix4.txt is the final matrix for FemXpress inference input that based on the SNP sites that are preserved under the strictest and most credible condition.

Step 2 : run FemXpress using the barcode-SNP matrix obtained in the previous step.For instance:
python FemXpress_2.py result_matrix4.txt > nohup_real.txt
