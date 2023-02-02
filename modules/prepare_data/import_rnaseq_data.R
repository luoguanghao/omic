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
library(DESeq2)

library(edgeR)
library(limma)

library(NetBID2)

# get VIZOME data #########################################
## mut data

maf_df = read.csv('/y/Archive/Bondi/data/VIZOME/LAML_nature_PDC.maf',sep='\t')
### 获取FLT3_ITD
ITD_sample_ls = as.vector((maf_df %>% filter(Hugo_Symbol=='FLT3'&ITDorNOT=='ITD'))['Tumor_Sample_Barcode'][['Tumor_Sample_Barcode']])
ITD_sample_ls = str_replace_all(paste('X',ITD_sample_ls,sep=''),'-','.') # colname的格式化



## exp data
exp_mat = read.csv('/y/Archive/Bondi/data/VIZOME/nature_aml_log2_fpkm.txt',sep='\t')
# 处理重复基因，滤出需要的基因
exp_mat = aggregate(.~ Symbol, exp_mat, mean)
exp_mat = exp_mat %>% column_to_rownames('Symbol')
sample_ls = colnames(exp_mat)


## drug sen data
drug_sen_df = read.csv('/y/Archive/Bondi/data/VIZOME/nature_aml_drug_sen.tsv',sep='\t')
drug_sen_df$lab_id = lapply(drug_sen_df$lab_id, function(x) str_replace_all(paste('X',x,sep=''),'-','.')) # colname的格式化
drug_sen_df$lab_id = as.factor(as.character(drug_sen_df$lab_id))
drug_ls = drug_sen_df[['inhibitor']][!duplicated(drug_sen_df[['inhibitor']])]


# get CH data #########################################


#########################################
# get TCGA data #########################
#########################################

## exp mat #########################################
exp_mat = read.csv('/y/Archive/Bondi/data/TCGA/LAML/TCGA-LAML.htseq_fpkm.tsv.gz',sep='\t')
probeMap_df = read.csv('/y/Archive/Bondi/data/TCGA/LAML/gencode.v22.annotation.gene.probeMap',sep='\t')
###
probeMap_df <- probeMap_df %>%  dplyr::select(1,2)
colnames(probeMap_df) <- c("Ensembl_ID","gene_name")

exp_mat <- inner_join(probeMap_df, exp_mat, by = "Ensembl_ID") %>% dplyr::select(-1)
exp_mat <- aggregate(.~ gene_name, exp_mat, mean)
# save()
exp_mat <- exp_mat %>% column_to_rownames("gene_name")

exp_mat = 2**(exp_mat)-1 # for tcga count data

## ========== new version

exp_mat = read_tsv(exp_path)
probeMap_df = read.csv(probeMap_path,sep='\t')
###
probeMap_df <- probeMap_df %>%  dplyr::select(1,2)
colnames(probeMap_df) <- c("Ensembl_ID","gene_name")

exp_mat[,c(2:dim(exp_mat)[2])] = 2**exp_mat[-1]-1 # for count data, reverse

exp_mat <- inner_join(probeMap_df, exp_mat, by = "Ensembl_ID") %>% dplyr::select(-1)
exp_mat <- aggregate(.~ gene_name, exp_mat, mean)
rn = rownames(exp_mat)
exp_mat <- exp_mat %>% column_to_rownames("gene_name")
exp_mat = apply(exp_mat,2,as.integer)
row.names(exp_mat) = rn

## mut #########################################

## methy #########################################


## survival #########################################



#########################################
# get depmap data #######################
#########################################

## sample info
info_df = read.csv('/y/Archive/Bondi/data/depmap/info/sample_info.csv',sep=',')



## exp mat
exp_df = read_csv('~/my_project/data/depmap/expression/CCLE_expression.csv')
colnames(exp_df) = sapply(strsplit(colnames(exp_df),' \\('),'[[',1)
colnames(exp_df)[1] = 'DepMap_ID'
exp_df = t(exp_df%>%column_to_rownames('DepMap_ID'))

## mut


## drug sen
drug_sen_df = read.csv('/y/Archive/Bondi/data/depmap/drug_sensitivity/sanger-dose-response.csv')
gdsc_cell_ls = drug_sen_df[['Cell.line.name']]

## metabolics


## get sp cell line
sp_cell_line = info_df[ grep('AML',info_df[['Subtype']]), ][['DepMap_ID']]
cb_exp_mat = merge(info_df['DepMap_ID'],exp_mat,by="DepMap_ID",all=F) %>% column_to_rownames('DepMap_ID') %>% t()
cb_exp_mat = cb_exp_mat[,intersect(sp_cell_line,colnames(cb_exp_mat))] # cell line不一定有转录组信息，选出有转录组信息的
















