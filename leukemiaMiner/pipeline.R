
source('/home/lgh/my_project/omic/modules/plot/survival.R')
source('/home/lgh/my_project/omic/modules/grouping/func_post_subtyping.R')
source('/home/lgh/my_project/omic/modules/feature_selection/by_var.R')

library(IOBR)
#library(EPIC)
#library(estimate) 
library(tidyverse)
library(tidyHeatmap)
library(maftools)
library(ggpubr)
library(ggplot2)
library(survival)
library(patchwork)

library(tidyverse)
library(DESeq2)
library(GSVA)
library(pheatmap)
library(ConsensusClusterPlus)
library(ggfortify)

# source('/home/lgh/my_project/omic/modules/prepare_data/depmap_data.R')
# source('/home/lgh/my_project/omic/project/protostatsis/func_subtyping.R')
source('/home/lgh/my_project/omic/modules/grouping/func_subtyping.R')
source('/home/lgh/my_project/omic/modules/bioinformatics_wrapper/gsva.R')
source('/home/lgh/my_project/omic/modules/grouping/grouping_by_clustering_pipeline.R')
source('/home/lgh/my_project/omic/modules/grouping/func_post_subtyping.R')

source('/home/lgh/my_project/omic/leukemiaMiner/discovery_cluster_analysis.R')
source('/home/lgh/my_project/omic/leukemiaMiner/discovery_factor_analysis.R')

source('/home/lgh/my_project/omic/leukemiaMiner/predict.R')
source('/home/lgh/my_project/omic/leukemiaMiner/utils.R')
library(clusterProfiler)


######
######

ch_exp_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_deseq2norm_expression_cleaned.tsv')
ch_surv_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_surv_rfs_eln_fab_0401_cleaned.tsv')
ch_exp_mat_scale = t(scale(t(ch_exp_df%>%column_to_rownames('Hugo_Symbol'))))

# ch_pro_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_tpm_expression_cleaned.tsv')
ch_drug_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_rawic50_drug_sensitivity_cleaned.tsv')
mut_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_mutations_matrix_cleaned.tsv')


vizome_exp_df = read_tsv('/home/lgh/my_project/data/VIZOME/cleaned/VIZOME_log2fpkm_expression_cleaned.tsv')
vizome_surv_df = read_tsv('/home/lgh/my_project/data/VIZOME/cleaned/VIZOME_survival_cleaned.tsv')
vizome_exp_mat_scale = t(scale(t(vizome_exp_df%>%column_to_rownames('Hugo_Symbol'))))


######## ch pro
# pho_mat = read_tsv('/home/lgh/my_project/data/changhai_data/protein_phosph/second/sort_by_sample/impute/new/phosphosite_df_min_impute_scalegene.tsv')%>%column_to_rownames('Hugo_Symbol')
pro_mat = read_tsv('/home/lgh/my_project/data/changhai_data/protein_phosph/second/sort_by_sample/impute/new/protein_df_min_impute_scalegene.tsv')%>%column_to_rownames('Hugo_Symbol')

pro_mat_t = t(pro_mat)

pho_mat_raw = read_tsv('/home/lgh/my_project/data/changhai_data/protein_phosph/second/sort_by_sample/CH_amlid_phosphoprotein_site_withna_nolog_geneXsample.tsv')%>%column_to_rownames('Hugo_Symbol')
pro_mat_raw = read_tsv('/home/lgh/my_project/data/changhai_data/protein_phosph/second/sort_by_sample/CH_amlid_protein_withna_nolog_geneXsample.tsv')%>%column_to_rownames('Hugo_Symbol')

sel_pro_by_na = rownames(pro_mat_raw)[apply(pro_mat_raw,1,function(x) sum(is.na(x)))<101*0.5]
sel_pho_by_na = rownames(pho_mat_raw)[apply(pho_mat_raw,1,function(x) sum(is.na(x)))<101*0.5]


######## AML2022 
pro_2022_mat = read_csv('/home/lgh/my_project/data/AML_Proteomics_2022/pro_impute_Discovery_Cohort.csv')
pro_2022_mat = cbind(pro_2022_mat[-1,179],pro_2022_mat[-1,-c(178:181)])

pro_2022_mat[,2:ncol(pro_2022_mat)] = sapply(pro_2022_mat[,2:ncol(pro_2022_mat)],as.numeric)
pro_2022_mat = aggregate(.~PG.Genes,pro_2022_mat,mean)
pro_2022_mat_t = data.frame(t(pro_2022_mat%>%column_to_rownames('PG.Genes')))
colnames(pro_2022_mat_t) = sapply(colnames(pro_2022_mat_t), function(x) if(length(strsplit(x,'\\.')[[1]])>1){strsplit(x,'\\.')[[1]][2]}else{x} )
pro_2022_mat_t = scale(pro_2022_mat_t)
                                  
pro2022_survival = read_tsv('/home/lgh/my_project/data/AML_Proteomics_2022/Clinical_Discovery_Cohort.txt')[,c('ID','Death Event','OS [months]','Relapse Event','RFS [months]')]
colnames(pro2022_survival) = c('sample','OSS','OS','RFSS','RFS')


pro_2022_valid_mat = read_csv('/home/lgh/my_project/data/AML_Proteomics_2022/pro_impute_Validation_Cohort.csv',skip=1)
pro_2022_valid_mat = cbind(pro_2022_valid_mat[-1,91],pro_2022_valid_mat[-1,-c(76:91)])
colnames(pro_2022_valid_mat)[1] = 'Gene names'
pro_2022_valid_mat[,2:ncol(pro_2022_valid_mat)] = sapply(pro_2022_valid_mat[,2:ncol(pro_2022_valid_mat)],as.numeric)
pro_2022_valid_mat = aggregate(.~`Gene names`,pro_2022_valid_mat,mean)
pro_2022_valid_mat_t = data.frame(t(pro_2022_valid_mat%>%column_to_rownames('Gene names')))
colnames(pro_2022_valid_mat_t) = sapply(colnames(pro_2022_valid_mat_t), function(x) if(length(strsplit(x,'\\.')[[1]])>1){strsplit(x,'\\.')[[1]][2]}else{x} )

pro2022_valid_survival = read_tsv('/home/lgh/my_project/data/AML_Proteomics_2022/Clinical_Validation_Cohort.txt')[,c('ID','Death Event','OS [months]','Relapse Event','RFS [months]')]
colnames(pro2022_valid_survival) = c('sample','OSS','OS','RFSS','RFS')
                                        
pro_2022_valid_mat_t = scale(pro_2022_valid_mat_t)


PATH_DATA = '/home/lgh/my_project/data'
gmtfile = list(
    gobp=paste(PATH_DATA,'/pathway/c5.go.bp.v7.5.1.symbols.gmt',sep=''),
    kegg=paste(PATH_DATA,'/pathway/c2.cp.kegg.v7.5.1.symbols.gmt',sep=''),
    reactome=paste(PATH_DATA,'/pathway/c2.cp.reactome.v7.5.1.symbols.gmt',sep=''),
    hallmark='/home/lgh/my_project/data/pathway/h.all.v7.5.1.symbols.gmt'
    #tft=paste(PATH_DATA,'/pathway/c3.tft.v7.5.1.symbols.add_cebpa.gmt',sep='')
)





























