library(stringr)
library(ggpubr) 
library(ggplot2)
library(GSEABase)
library(GSVAdata)
#data(c2BroadSets)
#c2BroadSets

library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)
library(tidyverse)
library(pheatmap)

library(ConsensusClusterPlus)


# pathway info #################
pathway_ls_path = '/y/Archive/Bondi/jupyter/hl/模仿Cell Meta分析/MetaPathway_ls.txt'
pw_data <- read.csv(pathway_ls_path, encoding="UTF-8",skip=0,sep='\t',stringsAsFactors = FALSE)

gene_vec = c()
pathway_vec = c()
for(i_r in 1:dim(pw_data)[1]){
# for(i_r in 1:3){
    gene = pw_data[i_r,1]
    pathways = as.character(unlist(strsplit(pw_data[i_r,2], split = "; ")))
    for(pw in pathways){
        gene_vec = c(gene_vec,gene)
        pathway_vec = c(pathway_vec,pw)
    }
}
gene_pathway_df = data.frame(gene=gene_vec, pathway=pathway_vec)
gene_pathway_list <- split(as.matrix(gene_pathway_df)[,1], gene_pathway_df[,2])

pw_g_list = split(as.matrix(gene_pathway_df)[,1], gene_pathway_df[,2])


# expression data #################
## from VIZOME
exp_mat = read.csv('/y/Archive/Bondi/data/VIZOME/nature_aml_log2_fpkm.txt',sep='\t')
exp_mat = aggregate(.~ Symbol, exp_mat, mean) #处理重复基因名问题
rownames(exp_mat) = exp_mat$Symbol
exp_mat = exp_mat[-1]

## from TCGA
exp_mat = read.csv('/y/Archive/Bondi/data/TCGA/LAML/TCGA-LAML.htseq_fpkm.tsv.gz',sep='\t')
probeMap_df = read.csv('/y/Archive/Bondi/data/TCGA/LAML/gencode.v22.annotation.gene.probeMap',sep='\t')

probeMap_df <- probeMap_df %>%  dplyr::select(1,2)
colnames(probeMap_df) <- c("Ensembl_ID","gene_name")

exp_mat <- inner_join(probeMap_df, exp_mat, by = "Ensembl_ID") %>% dplyr::select(-1)
exp_mat <- aggregate(.~ gene_name, exp_mat, mean)
exp_mat <- exp_mat %>% column_to_rownames("gene_name")

# GSVA #################
# exp_mat: gene×sample
##
exp_mat = exp_mat[pw_data$Gene.symbol,] # 选择pathway基因
exp_mat = na.omit(exp_mat)

gsva_matrix <- gsva(as.matrix(exp_mat), pw_g_list,method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE, parallel.sz=1)

# Consensus Cluster #################
# input of ConsensusClusterPlus: gene×sample ,this func is do the cluster for sample
##
results = ConsensusClusterPlus(t(scale(t(gsva_matrix))),maxK=8)

save(results, file = "./result/consenCluster_res_vizome.RData")


## post cluster analysis

pw_sort_ls = read.table('./pathway_sort_ls.txt',header = FALSE,sep = '\n')

### 这里是取4类的结果来分析
class_df = as.data.frame(results[[4]][['consensusClass']],col.names=c('class'))
colnames(class_df)=c('class')
class_df['sample'] = rownames(class_df)

tmp1 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['1']]
tmp2 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['2']]
tmp3 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['3']]
tmp4 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['4']]

### 取6类的结果来分析
class_df = as.data.frame(results[[6]][['consensusClass']],col.names=c('class'))
colnames(class_df)=c('class')
class_df['sample'] = rownames(class_df)

tmp1 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['1']]
tmp2 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['2']]
tmp3 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['3']]
tmp4 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['4']]
tmp5 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['5']]
tmp6 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['6']]

### 画热图

#### 4c
annotation_col = data.frame(
  Cluster = factor(c(rep('c1',length(tmp1)) , rep('c2',length(tmp2)) , rep('c3',length(tmp3)) , rep('c4',length(tmp4))))
)
rownames(annotation_col) = c(tmp1,tmp2,tmp3,tmp4)

rownames(gsva_matrix) = toupper(as.character(rownames(gsva_matrix)))
pheatmap(
    gsva_matrix[toupper(as.character(pw_sort_ls[[1]])),c(tmp1,tmp2,tmp3,tmp4)],
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    scale = 'row',
    color = colorRampPalette(c("green","black", "red"))(50),
    #filename='./sorted_pw_cencluster_hmap3.pdf',
    cellwidth = 2, cellheight = 15,
    annotation_col = annotation_col,
    border_color = "white"
    )

#### 6c
pw_sort_ls = read.table('./pathway_sort_ls.txt',header = FALSE,sep = '\n')
rownames(gsva_matrix) = toupper(as.character(rownames(gsva_matrix)))

annotation_col = data.frame(
  Cluster = factor(c(rep('c1',length(tmp1)) , rep('c2',length(tmp2)) , 
                     rep('c3',length(tmp3)) , rep('c4',length(tmp4)),
                     rep('c5',length(tmp5)) , rep('c6',length(tmp6)))
                )
)
rownames(annotation_col) = c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)

pheatmap(
    gsva_matrix[toupper(as.character(pw_sort_ls[[1]])),c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)],
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    scale = 'row',
    color = colorRampPalette(c("green","black", "red"))(50),
    filename='./result/sorted_pw_cencluster_hmap_vizome_6c.pdf',  ## 改名
    cellwidth = 3, cellheight = 20,
    annotation_col = annotation_col,
    border = FALSE
    # border_color = "white"
    )

































