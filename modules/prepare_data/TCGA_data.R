library(tidyverse)
library(stringr)

#############
# rnaseq
#############
probeMap_dir = '/y/Archive/Bondi/data/TCGA/LAML/gencode.v22.annotation.gene.probeMap'
exp_dir = '/y/Archive/Bondi/data/VIZOME/nature_aml_log2_fpkm.txt'
###
probeMap_df = read.csv(probeMap_dir,sep='\t')
exp_mat = read.csv(exp_dir,sep='\t')
# probeMap 标注


probeMap_df <- probeMap_df %>%  dplyr::select(1,2)
colnames(probeMap_df) <- c("Ensembl_ID","gene_name")
exp_mat <- inner_join(probeMap_df, exp_mat, by = "Ensembl_ID") %>% dplyr::select(-1)
# 处理重复基因，滤出需要的基因
exp_mat <- aggregate(.~ gene_name, exp_mat, mean)
exp_mat <- exp_mat %>% column_to_rownames("gene_name")

###
#exp_mat
#############
#############



#############
# mut
#############


#############
#############



#############
# sample_data
#############


#############
#############

#############
# survival
#############
read_tsv('/home/lgh/my_project/data/TCGA/LAML/TCGA-LAML.survival.tsv.gz')

#############
#############

