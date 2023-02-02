# http://localhost:9424/notebooks/my_home/my_project/jn/%E5%8F%8C%E7%AA%81%E5%8F%98%E4%B8%8B%E8%8D%AF%E6%95%8F%E6%AF%94%E8%BE%83.ipynb
#


library(tidyverse)
library(ggpubr)

drug_sen_df = read_csv('~/my_project/data/CCLE/drug_sensitivity/CCLE_NP24.2009_Drug_data_2015.02.24.csv')
maf_df = read_tsv('~/my_project/data/CCLE/genome_data/CCLE_DepMap_18q3_maf_20180718.txt')
# 药物信息，包含靶点信息
drug_info = read_csv('~/my_project/data/depmap/info/screened_compunds_rel_8.2.csv')

# get pi3k drug
pi3k_ls = toupper((drug_info%>%filter(grepl('PI3K',TARGET)))$DRUG_NAME)

# 建立df描述每个样本突变情况
tmp_mut_df = data.frame(DepMap_ID = unique(maf_df$DepMap_ID))
tmp_mut_df[['TP53']] = 0
tmp_mut_df[['PIK3CA']] = 0

mut_sample_ls = (maf_df%>%filter(Hugo_Symbol=='TP53'))$DepMap_ID
tmp_mut_df[tmp_mut_df$DepMap_ID%in%mut_sample_ls,]$TP53 = 1
mut_sample_ls = (maf_df%>%filter(Hugo_Symbol=='PIK3CA'))$DepMap_ID
tmp_mut_df[tmp_mut_df$DepMap_ID%in%mut_sample_ls,]$PIK3CA = 1

# 加入药敏信息
tmp_drugsen_mut_df = merge(drug_sen_df[,c('ARXSPAN_ID','auc','DRUG_NAME')], 
                            tmp_mut_df, by.x='ARXSPAN_ID', by.y='DepMap_ID')
tmp_drugsen_mut_df = na.omit(tmp_drugsen_mut_df%>%filter(DRUG_NAME%in%pi3k_ls))

# 加入描述共突变情况
tmp_drugsen_mut_df$mut_type='wt'
tmp_drugsen_mut_df[tmp_drugsen_mut_df$TP53==1,]$mut_type='TP53'
tmp_drugsen_mut_df[tmp_drugsen_mut_df$PIK3CA==1,]$mut_type='PIK3CA'
tmp_drugsen_mut_df[tmp_drugsen_mut_df$TP53==1&tmp_drugsen_mut_df$PIK3CA==1,]$mut_type='both'

# plot
comp = as.data.frame(combn(unique(tmp_drugsen_mut_df$mut_type),2))
my_comparisons = list(comp$V1, comp$V2, comp$V3, comp$V4, comp$V5, comp$V6)

p <- ggboxplot(tmp_drugsen_mut_df, 
          x = "mut_type", y = "auc",
          color = "mut_type", palette = "nature", facet.by = "DRUG_NAME")+
          theme(axis.text.x = element_text(angle = 75, vjust=0.5)) +
          stat_compare_means(aes(group=mut_type),comparisons=my_comparisons)
          #add = "jitter")
ggsave(sprintf('out2.pdf'), dpi=50, width=15, height=15)


