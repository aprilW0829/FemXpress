# FemXpress
FemXpress--tool that classifies the source of X chromosome inactivation of each female cell.

## User manual and guide
### Overview
FemXpress is a Python package, a powerful tool designed for a comprehensive analysis of XCI status in female scRNA-Seq datasets. It has the following main functions: 
1) to infer the source of X chromosome inactivation (paternal or maternal) for each cell,
2) to infer new potential escape genes,
3) to infer the X-to-auto expression dose in different cell types.

#### Workflow of FemXpress
To avoid relying on parental SNP information, we performed de novo discovery of SNP that exhibit parental origin differences using the alignment results of the pseudo-bulk scRNA-Seq sample. In order to minimize the impact of sequencing and analysis errors on SNP calling, we applied strict filtering criteria, including: 
1) a minimum coverage of at least 2 reads supporting either of the two genotypes of the SNP,
2) a minimum number of 2 cells with each genotype,
3) exclusion of loci originating from simple repeat regions,
4) sorting the remaining potential SNPs based on the ratio of reads for each genotype at the locus and removing the top and bottom 10% of SNPs with the highest and lowest ratios, respectively, as they are likely introduced by sequencing errors in high coverage regions (Fig a-c).

Subsequently, linkage analysis was performed based on the identified high-confidence SNPs. We counted the genotypes of two SNPs that appear simultaneously in each cell, considering them linked genotypes only if both genotypes are present in at least 80% of cells (Fig d). Next, we conducted an extended linkage scoring by combining linked genotypes into longer chains of linked loci. The pair of haplotypes with the highest number of linked loci was selected as the genotypes of the parental source X chromosome (Fig e). Then, when the haplotypes of the parents were phased, the origin of the X chromosome for each cell was determined by a voting method. Each haplotype received one vote if the cell carried an allele base from that haplotype, and the cell was classified to have that haplotype if the votes were higher than those of the other haplotype. Finally, in the cells that were already classified with the same parental source of the inactivated X chromosome, we examined the frequency of concurrent parental origin genotypes at each SNP locus to identify potential XCI escaping gene derived SNPs(Fig f).
![image](https://github.com/wangxin970829/FemXpress/blob/main/images/workflow.jpg)


## Citing FemXpress
A preprint FemXpress: Systematic Analysis of X Chromosome Inactivation Heterogeneity in Female Single-Cell RNA-Seq Samples provides a detailed overview of the statistical models used in FemXpress. We ask that you cite this paper if you use FemXpress in work that leads to publication.This preprint is used for documentation and citation.

## Install
To install FemXpress, make sure you need install some dependances, if you need more details on the dependences, look at the FemXpress-env.yaml file.
set up conda environment for FemXpress
```
conda env create -f FemXpress-env.yaml
```

Install FemXpress from zip package
```
conda activate FemXpress-env
wget https://github.com/wangxin970829/FemXpress/archive/refs/heads/main.zip
unzip main.zip
cd FemXpress-main/scripts
export PATH=$PATH:/path/to/your/FemXpress-main/scripts/FemXpress_1.py
export PATH=$PATH:/path/to/your/FemXpress-main/scripts/FemXpress_2.py
```

## Usage
The current version of FemXpress was only developed based on single-cell sequencing data from 10×Genomics. It consists of two running steps:

Step 1 : to obtain pre-processed barcode-SNP matrix based on the alignment bam file from cellranger’s output.
According to how it works, it needs to have the following input files: BAM file for CellRanger alignment output, cell tag file (which can be meta information file that has been processed by Seurat)(possorted_genome_bam.bam), reference genome sequence file(.fa), and rmsk file.For instance:

Run the following arguments for command-line help:
```
$ python FemXpress_1.py --help
usage: FemXpress_1.py [-h] [-b BAM] [-g GENOME] [-m META] [-r RMSK] [-v]

preprocess of FemXpress, Example: python FemXpress_1.py -b possorted_genome_bam.bam -g genome.fa -e meta.txt -r rmsk.txt

optional arguments:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     input bam file from CellRanger output
  -g GENOME, --genome GENOME
                        input genome
  -m META, --meta META  input meta information of barcodes
  -r RMSK, --rmsk RMSK  rmsk file
  -v, --verbose         print verbose output
```

...will produce four barcode-SNP matrix files in the subdirectory ./FemXpress/result, result_matrix4.csv is the final matrix for FemXpress inference input that based on the SNP sites that are preserved under the strictest and most credible condition.

Step 2 : run FemXpress using the barcode-SNP matrix obtained in the previous step.For instance:
```
$ python FemXpress_2.py --help
usage: FemXpress_2.py [-h] [-m MATRIX] [-a ANNOT] [-v]

inference of FemXpress, Example: python FemXpress_2.py -m result_matrix4.csv -a /data2/wangxin/database/genome/gencode.vM25.chr_patch_hapl_scaff.basic.annotation.gtf >
nohup_real.txt

optional arguments:
  -h, --help            show this help message and exit
  -m MATRIX, --matrix MATRIX
                        output result_matrix4.csv generated in the previous step
  -a ANNOT, --annot ANNOT
                        gene annotation txt[gtf]
  -v, --verbose         print verbose output

```

## Output
The output will generate three subdirectories in the ./FemXpress directory, which contains tmp, result, inference.  

The tmp subdirectory contains the generated intermediate files such as histograms of the statistical distribution of the SNPs obtained at the 4 thresholds,

The result subdirectory generates four matrix files(barcodes in rows,SNP sites in columns) corresponding to the four SNP sites file under the ./FexmXpress/tmp, for example:
|                  |chrX;159974518|chrX;147053550|chrX;17564774|chrX;47919892|......|
|:-:|:-:|:-:|:-:|:-:|:-:|
|CCTTCAGGTGTAAATG-1|.;.;.;.|.;.;.;.|.;.;.;.|4;0;0;1|......|
|TTCCTAAGTCGTCATA-1|2;0;0;0|.;.;.;.|.;.;.;.|.;.;.;.|......|
|CCTTCAGGTGTAAATG-1|.;.;.;.|.;.;.;.|.;.;.;.|.;.;.;.|......|
|CCTTCAGGAGTACATG-1|.;.;.;.|0;3;1;0|.;.;.;.|.;.;.;.|......|
|......|......|......|......|......|......|

The inference subdirectory is generated as an inferred result for each cell by FemXpress, it contains 18 files, of which escape_pos_list_method_3.txt is the inferred new potential escape genes, clusters_vote_method_1.tsv is the classification result inferred(0 and 1 are the main two categories, and 2 is not successfully classified by FemXpress). For example:
escape_pos_list_method_3.txt
|151869594|
|:-:|
|159646980|
|7307213|
|105052403|
|10434337|

clusters_vote_method_1.tsv
|GTTCGCTGTACGGCAA-1|0|
|:-:|:-:|
|GACTCTCAGTGGATAT-1|0|
|TCCAGAAGTCGTGATT-1|1|
|CAATCGACAGTTCCAA-1|0|
|CAAGCTATCTGCTTAT-1|0|




