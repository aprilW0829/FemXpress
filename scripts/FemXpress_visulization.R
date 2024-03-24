#! /usr/lib/R/bin/Rscript --vanilla

library(optparse)
option_list <- list(#构建参数列表
    make_option( c('-i', '--sample'), type = 'character', default = FALSE, help = 'sample name'),  
    make_option( c('-r', '--rds'), type = 'character', default = FALSE, help = 'seurat_rds'),
    make_option( c('-c', '--cluster'), type = 'character',default = FALSE, help = 'cluster.tsv'),
    make_option( c('-x', '--chrX_genes'),type = 'character',default = FALSE,help = 'chrX_genes_txt')
    )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
sample1 = opt$sample
rds1 = opt$rds
cluster1 = opt$cluster
chrX_genes1 = opt$chrX_genes

visulizaation_dir = './FemXpress/visulization'

Afunction <- function(sample,rds,cluster,chrX_genes){
print('load needed packages!')
library(optparse)
library(Seurat)
library(SingleR)
library(utils)
library(stringr)
library(dplyr)
library(tidyverse)
library(patchwork)
library(harmony)
library(clustree)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(SummarizedExperiment)
library(scuttle)
library(reshape2)
library(readxl)
library(ggpubr)

print(str_c(sample, rds, cluster,chrX_genes))

#human_data <- read_excel('/data5/wangxin/20220923fgh_mouse_brain/puplic_data/human_data_info/subclusterX_in_human_application.xlsx',sheet=8)
#human_data <- as.data.frame(human_data)
#human_data_gsm <- human_data[!duplicated(human_data$GSM),]
#gsm_sample <- read.table('/data5/wangxin/20220923fgh_mouse_brain/puplic_data/subclusterX/GSM_sample.txt',sep='\t',header=TRUE)

#print(str_c(rds))
rds <- readRDS(rds)
#print(rds)
cluster <- read.table(cluster,sep='\t',header=FALSE)
#print(cluster)
rds_sub <- subset(rds,cells=cluster$V1)

rds_sub@meta.data$tool_group <- cluster$V2
#rds_sub@meta.data$state_group <- subset(human_data_gsm,human_data_gsm$GSM==sample)$class

#class=subset(human_data_gsm,human_data_gsm$GSM==sample)$class

#saveRDS(rds_sub,file=str_c('/data5/wangxin/20220923fgh_mouse_brain/puplic_data/rds/',sample,'_tool.rds'))
anno_tool_group_num <- as.data.frame(table(rds_sub$anno,rds_sub$tool_group))
anno_tool_group_num_01 <- subset(anno_tool_group_num,anno_tool_group_num$Var2=='0'|anno_tool_group_num$Var2=='1')
#anno_tool_group_num_01$sample <- sample

# subclusterX 给每个样本中的细胞类型分类(0/1)的占比统计
p1 <- ggplot(data=anno_tool_group_num_01,aes(Var1,Freq,fill=Var2))+geom_bar(stat="identity",position="fill", color="black", width=0.7,size=0.25)+theme_bw()+theme_classic()+theme_bw()+theme(panel.grid.major=element_line(colour=NA))+theme(panel.grid.minor=element_line(colour=NA))+theme_classic()+theme(axis.title.x = element_blank(),axis.title.y=element_text(size = 10),axis.text.y = element_text(size = 10),axis.text.x = element_text(size = 10,angle=45, hjust=1, vjust=1))+scale_y_continuous(expand=c(0,0))+ylab('% of cells')+theme(legend.title=element_blank(),legend.text=element_text(size = 10))

ggsave(p1,file=str_c(visulization_dir,'/',sample,'_subclusterX_',class,'_barplot.pdf'),width=8,height=6)
#print(str_c(sample,' done!'))
write.table(anno_tool_group_num_01,file=str_c(visulization_dir,'/',sample,'_subclusterX_',class,'_subclusterX.txt'),sep='\t')


# 每个样本中不同细胞类型的占比统计
sample_num <- table(rds$orig.ident)
sample_anno_num <- table(rds$anno)
sample_anno_ratio <- as.data.frame(sample_anno_num/as.numeric(sample_num))
# celltypecolor <- c('F8766D','E68613','CD9600','ABA300','7CAE00','0CB702','00BE67','00C19A','00BFC4','00B8E7','00A9FF','8494FF','C77CFF','ED68ED','FF61CC')
pdf(str_c(visulization_dir,'/',sample,'_subclusterX_',class,'_anno_ratio_pie.pdf'),width=6,height=6)
pie(sample_anno_ratio$Freq, labels = sample_anno_ratio$Var1,col = rainbow(length(sample_anno_ratio$Freq)))
legend("right",sample_anno_ratio$Var1, cex = 0.8, fill = rainbow(length(sample_anno_ratio$Freq)),legend = "Label")
dev.off()
print(str_c(sample,' done!'))


# 每个样本中被推断出来的细胞占比统计
dev.new()
a <- table(rds_sub@meta.data$orig.ident)
b <- table(rds@meta.data$orig.ident)
infer_ratio <- as.numeric(a)/as.numeric(b)
uninfer_ratio <- 1-infer_ratio
ratio_data <- c(infer_ratio,uninfer_ratio)
labels <- c('inferred','uninferred')
pdf(str_c(visulization_dir,'/',sample,'_subclusterX_',class,'_infer_uninfer_ratio_pie.pdf'),width=6,height=6)
pie(ratio_data, labels = labels,col=c('F8766D','00BFC4'))
legend("right",sample_anno_ratio$Var1, cex = 0.8, fill = rainbow(length(sample_anno_ratio$Freq)))
dev.off()
print(str_c(sample,' done!'))


# 每个样本的细胞注释umap图
 dev.new()
 p1 <- DimPlot(rds,group.by='anno')
 ggsave(p1 ,file=str_c(visulization_dir,'/',sample,'_',class,'_singleR_anno_umap.pdf'),width=7,height=6)
 print('get anno umap for ',sample,' !')
 dev.off()


# get ratio(expression of genes that located on the chrX/genes that located on the autosomal chromosomes) of different celltyeps
a <- readRDS(rds)
a <- rds
chrx_gene <- read.table(chrX_genes,sep='\t',header=FALSE)
chrx_count <- a@assays$RNA@counts[chrx_gene$V1[chrx_gene$V1 %in% rownames(a@assays$RNA@counts)],]
chrx_count_t <- as.data.frame(t(chrx_count))
chrx_count_t$chrx_count_sum <- rowSums(chrx_count_t)

a$chrx_count_sum <- chrx_count_t$chrx_count_sum
a$chrauto_count <- a$nCount_RNA - a$chrx_count_sum
a$chrx_chrauto_count_ratio <- a$chrx_count_sum/a$chrauto_count

data <- a@meta.data[,c('chrx_chrauto_count_ratio','anno')]

counts <- table(data$anno)
result <- data[data$anno %in% names(counts[counts > 1]), ]
#print(result)

sample_id_name <- sample
p2 <- ggplot(result, aes(x=anno, y=chrx_chrauto_count_ratio, fill=anno)) + geom_boxplot() + theme_bw()+stat_compare_means(method = "t.test",label = "p.signif")+stat_boxplot(color='black',geom="errorbar")+theme_classic()+theme(panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank(),axis.title.x=element_text())+theme(legend.text=element_text(),legend.title=element_text())+xlab('go term')+ylab('count ratio of chrX/chrauto')+theme(axis.title.x = element_text(size=10),axis.title.y=element_text(size = 10),axis.text.y = element_text(size = 10),axis.text.x = element_text(size = 10,angle=45, hjust=1, vjust=1))+scale_y_continuous(expand=c(0,0))

ggsave(p2,file=str_c(visulization_dir,'/',sample_id_name,'_chrx_chrauto_count_ratio_barplot.pdf'),width=10,height=6)


write.table(result,file=str_c(visulization_dir,'/',sample_id_name,'_chrx_chrauto_ratio.txt'),sep='\t')
}

##### 调用函数并将运行结果保存输出
Afunction(sample1,rds1,cluster1,chrX_genes1)
