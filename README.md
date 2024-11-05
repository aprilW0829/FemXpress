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
![image](https://github.com/wangxin970829/FemXpress/blob/main/images/fig1.jpg)


## Install
You need to create a new environment for FemXpress
```
conda create -n FemXpress-env python=3.8
```
!Note We have put the binary compatibility tools including [freebayes](https://github.com/freebayes/freebayes), [samtools](https://github.com/samtools/samtools), [bam-readcount](https://github.com/genome/bam-readcount) in the apps folder. We fixed the version because the output formats vary a lot with different versions. If you are not able to run them, you can compile them in you system. We only test on these tools on following versions:

samtools version: 1.20 (using htslib 1.20)

bam-readcount version: v1.0.2

freebayes version: v1.0.2/v1.3.6


And in addition to these packages above, you also need some python dependences, you can refer to the documentation for details of these packages[requirements.txt] (https://github.com/wangxin970829/FemXpress/blob/main/requirements.txt) following the steps below.

Set up conda environment and install the required python packages for FemXpress
```
wget https://github.com/wangxin970829/FemXpress/archive/refs/heads/main.zip
unzip main.zip
cd FemXpress-main/
conda activate FemXpress-env
pip install -r requirements.txt
chmod 755 apps/*
cd scripts
path="XXX/FemXpress-main"  # where FemXpress is downloaded
python ${path}/scripts/FemXpress_1.py --help
python ${path}/scripts/FemXpress_2.py --help
```


## Usage
The current version of FemXpress was only developed based on single-cell sequencing data from 10×Genomics. The details of the input of FemXpress can be found in the [input](https://github.com/wangxin970829/FemXpress/tree/main/test/input)

# You can run FemXpress using test data as follows:
FemXpress_PATH={path_to_FemXpress-main}
cd ${FemXpress_PATH}/test/
unzip *gz
python ../script/FemXpress_1.py -b ./input/test200000.bam -g ./input/mm10_chrX.fa -e ./input/meta_half.txt -r ./input/mm10_rmsk_test.txt -n test8.0 -c 2 > nohup_real1.txt
python ../script/FemXpress_2.py -m ./FemXpress/test8.0/result/result_matrix1.csv -a ./input/gencode.v44.chr_patch_hapl_scaff.basic.annotation.gtf -n test8.0 > nohup_real2.txt


It consists of two running steps:
### Step 1: 
To obtain pre-processed barcode-SNP matrix based on the alignment bam file from cellranger’s output.
According to how it works, it needs to have the following input files: BAM file for CellRanger alignment output(possorted_genome_bam.bam), cell tag file (which can be meta information file that has been processed by Seurat, barcodes's format located in the first column need to be the format which is the same as the bam' CB tag, such as:ACGCGCGCTACGCGCT-1), reference genome sequence file(.fa), and rmsk file(the five and six columns need to be the location of rmsk sequencess). For instance:

Run the following arguments for command-line help:
```
$ python ${path}/scripts/FemXpress_1.py --help
usage: FemXpress_1.py [-h] [-b BAM] [-g GENOME] [-m META] [-r RMSK] [-n sample] [-v]

preprocess of FemXpress, Example: python FemXpress_1.py -b possorted_genome_bam.bam -g genome.fa -e meta.txt -r rmsk.txt -n sample -c 2

optional arguments:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     input bam file from CellRanger output
  -g GENOME, --genome GENOME
                        input genome
  -e META, --meta META  input meta information of barcodes
  -r RMSK, --rmsk RMSK  rmsk file
  -n SAMPLE, --sample SAMPLE
                        sample name
  -c MINIMUM_COUNTS, --minimum_counts MINIMUM_COUNTS
                        minimum counts supporting second largest base, default: 2
  -v, --verbose         print verbose output
```

...will produce four barcode-SNP matrix files in the subdirectory ./FemXpress/result, result_matrix4.csv is the final matrix for FemXpress inference as input that based on the SNP sites that are preserved under the strictest and most credible condition.

### Step 2: 
Run FemXpress using the barcode-SNP matrix obtained in the previous step. For instance:
```
$ python ${path}/scripts/FemXpress_2.py --help
usage: python FemXpress_2.py [-h] [-m MATRIX] [-a ANNOT] [-n sample] [-v]

inference of FemXpress, Example: python FemXpress_2.py -m result_matrix4.csv -a /data2/wangxin/database/genome/gencode.vM25.chr_patch_hapl_scaff.basic.annotation.gtf -n sample > nohup_real.txt

optional arguments:
  -h, --help            show this help message and exit
  -m MATRIX, --matrix MATRIX
                        output result_matrix4.csv generated in the previous step
  -a ANNOT, --annot ANNOT
                        gene annotation txt[gtf]
  -n SAMPLE, --sample SAMPLE
                        sample name
  -v, --verbose         print verbose output
```

...will produce multiple files in the subdirectory ./FemXpress/inference, clusters_vote_method_1.tsv is the final inference file for each cell, escape_pos_list_method_3.txt is the inferenced new potential escape sites, the genes in which these sites are located are potential and new escape genes inferenced. And there are several other tmp files, for example: genotypes.tsv is the sites that are used to establish the chain relationship.
 

It provides one visulization step:
You can visulization FemXpress's inference result using [the procedures](https://github.com/wangxin970829/FemXpress/blob/main/scripts/FemXpress_visulization.R) and files produced in the previous step. For instance:
```
$ Rscript ${path}/scripts/FemXpress_visulization.R --help
Usage: subclusterX_visulization.R [options]

Options:
        -i SAMPLE, --sample=SAMPLE
                sample name
        -r RDS, --rds=RDS
                seurat_rds
        -c CLUSTER, --cluster=CLUSTER
                cluster.tsv
        -x CHRX_GENES, --chrX_genes=CHRX_GENES
                chrX_genes_txt
        -h, --help
                Show this help message and exit

Example: FemXpress_visulization.R -i sample -r rds -c clusters.tsv -x chrX_genes.txt
```
Or you can choose to write codes to make more beautiful pdf.


## Output
The output will generate three subdirectories in the ./FemXpress directory, which contains tmp, result, inference，visulization.  The details of the ouput of FemXpress can be found in the [output](https://github.com/wangxin970829/FemXpress/tree/main/test/output)

The [tmp](https://github.com/wangxin970829/FemXpress/tree/main/test/output/tmp) subdirectory contains the generated intermediate files such as histograms of the statistical distribution of the SNPs obtained at the 4 thresholds,

The [result](https://github.com/wangxin970829/FemXpress/tree/main/test/output/result) subdirectory generates four matrix files(barcodes in rows,SNP sites in columns) corresponding to the four SNP sites file under the ./FexmXpress/tmp, for example:

result_matrix4.csv
|                  |chrX;159974518|chrX;147053550|chrX;17564774|chrX;47919892|......|
|:-:|:-:|:-:|:-:|:-:|:-:|
|CCTTCAGGTGTAAATG-1|.;.;.;.|.;.;.;.|.;.;.;.|4;0;0;1|......|
|TTCCTAAGTCGTCATA-1|2;0;0;0|.;.;.;.|.;.;.;.|.;.;.;.|......|
|CCTTCAGGTGTAAATG-1|.;.;.;.|.;.;.;.|.;.;.;.|.;.;.;.|......|
|CCTTCAGGAGTACATG-1|.;.;.;.|0;3;1;0|.;.;.;.|.;.;.;.|......|
|......|......|......|......|......|......|

The [inference](https://github.com/wangxin970829/FemXpress/tree/main/test/output/inference) subdirectory is generated as an inferred result for each cell by FemXpress, it contains 18 files, of which escape_pos_list_method_3.txt is the inferred new potential escape genes, clusters_vote_method_1.tsv is the classification result inferred(0 and 1 are the main two categories, and 2 is not successfully classified by FemXpress). For example:

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

The [visulization](https://github.com/wangxin970829/FemXpress/tree/main/test/output/visulization) subdirectory generates four pdf files and two text files used for visualization.




