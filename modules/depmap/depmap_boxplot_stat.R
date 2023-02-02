# 这里是看在TP53 damaging背景下，RB1是否damaging对USP7 ceres score的影响

library(tidyverse)
library(ggpubr)


crispr_df = read_csv('~/my_project/data/depmap/dependence/CRISPR_gene_effect.csv')
colnames(crispr_df) = sapply(strsplit(colnames(crispr_df),' \\('),'[[',1)
mut_df = read_csv('~/my_project/data/depmap/genome/CCLE_mutations.csv')


df_for_plot = crispr_df[,c('DepMap_ID','USP7')]
df_for_plot$MDM2 = '0'
df_for_plot$TP53 = '0'
df_for_plot$RB1 = '0'


df_for_plot[df_for_plot$DepMap_ID%in%(mut_df%>%filter(Hugo_Symbol=='TP53'&Variant_annotation=='damaging'))$DepMap_ID,]$TP53 = 'damaging'
df_for_plot[df_for_plot$DepMap_ID%in%(mut_df%>%filter(Hugo_Symbol=='RB1'&Variant_annotation=='damaging'))$DepMap_ID,]$RB1 = 'damaging'
df_for_plot[df_for_plot$DepMap_ID%in%(mut_df%>%filter(Hugo_Symbol=='MDM2'&Variant_annotation=='damaging'))$DepMap_ID,]$MDM2 = 'damaging'


ggviolin(df_for_plot%>%filter(TP53=='damaging'),
        x='RB1',y='USP7',color='RB1',add='jitter')+
        stat_compare_means()

ggstatsplot::ggbetweenstats(df_for_plot%>%filter(TP53=='damaging'),
         x='RB1',y='USP7',color='RB1',add='jitter',type='nonparametric')




















