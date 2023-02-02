library(data.table)
library(tidyverse)
library(ggplot2)

library(ggpubr)
library(ggbeeswarm)
library(ggthemes)

options (warn = - 1)

## mut data
maf_df = as.data.frame(fread('~/my_project/data/VIZOME/LAML_nature_PDC.maf',sep='\t'))
maf_df$Tumor_Sample_Barcode = str_replace_all(paste('X',maf_df$Tumor_Sample_Barcode,sep=''),'-','.')

ITD_sample_ls = as.vector((maf_df %>% dplyr::filter(Hugo_Symbol=='FLT3'&ITDorNOT=='ITD'))['Tumor_Sample_Barcode'][['Tumor_Sample_Barcode']])

exp_mat = read.table('~/my_project/data/VIZOME/nature_aml_log2_fpkm_cleaned.txt', sep='\t');head(exp_mat)


## drug sen data
drug_sen_df = as.data.frame(fread('~/my_project/data/VIZOME/nature_aml_drug_sen.tsv',sep='\t'))
drug_sen_df$lab_id = sapply(drug_sen_df$lab_id, function(x) str_replace_all(paste('X',x,sep=''),'-','.'))
drug_sen_df$lab_id = as.factor(as.character(drug_sen_df$lab_id))
drug_ls = drug_sen_df[['inhibitor']][!duplicated(drug_sen_df[['inhibitor']])]


merge(maf_df%>%filter(Tumor_Sample_Barcode%in%ITD_sample_ls),
    drug_sen_df%>%filter(lab_id%in%ITD_sample_ls),by.x='Tumor_Sample_Barcode', by.y='lab_id')











