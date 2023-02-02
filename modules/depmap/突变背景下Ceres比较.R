library(tidyverse)
library(ggpubr)


exp_df = read_csv('~/my_project/data/depmap/expression/CCLE_expression.csv')
colnames(exp_df)[1] = 'DepMap_ID'
colnames(exp_df) = sapply(strsplit(colnames(exp_df),' \\('),'[[',1)

crispr_df = read_csv('~/my_project/data/depmap/dependence/CRISPR_gene_effect.csv')
colnames(crispr_df) = sapply(strsplit(colnames(crispr_df),' \\('),'[[',1)

mut_df = read_csv('~/my_project/data/depmap/genome/CCLE_mutations.csv')



# 比较突变背景下突变影响ceres
# example notebook：
df_for_plot = crispr_df[,c('DepMap_ID','USP7')]
df_for_plot$MDM2 = '0'
df_for_plot$TP53 = '0'
df_for_plot$RB1 = '0'

df_for_plot[df_for_plot$DepMap_ID%in%(mut_df%>%filter(Hugo_Symbol=='TP53'&Variant_annotation=='damaging'))$DepMap_ID,]$TP53 = 'damaging'
df_for_plot[df_for_plot$DepMap_ID%in%(mut_df%>%filter(Hugo_Symbol=='RB1'&Variant_annotation=='damaging'))$DepMap_ID,]$RB1 = 'damaging'
df_for_plot[df_for_plot$DepMap_ID%in%(mut_df%>%filter(Hugo_Symbol=='MDM2'&Variant_annotation=='damaging'))$DepMap_ID,]$MDM2 = 'damaging'

ggviolin(df_for_plot%>%filter(TP53=='damaging'),
         x='RB1',y='USP7',color='RB1',add='jitter')+stat_compare_means()
ggviolin(df_for_plot%>%filter(TP53!='damaging'),
         x='RB1',y='USP7',color='RB1',add='jitter',title='tp53 no damaging')+stat_compare_means()

ggstatsplot::ggbetweenstats(df_for_plot%>%filter(TP53=='damaging'),
         x='RB1',y='USP7',color='RB1',add='jitter',type='nonparametric')
ggstatsplot::ggbetweenstats(df_for_plot%>%filter(TP53!='damaging'),
         x='RB1',y='USP7',color='RB1',add='jitter',type='nonparametric')



# 比较突变背景下exp与ceres相关性
# example notebook：
df_for_plot = crispr_df[,c('DepMap_ID','USP7')]
#df_for_plot$MDM2 = '0'
df_for_plot$TP53 = '0'
#df_for_plot$RB1 = '0'

df_for_plot[df_for_plot$DepMap_ID%in%(mut_df%>%filter(Hugo_Symbol=='TP53'&Variant_annotation=='damaging'))$DepMap_ID,]$TP53 = 'damaging'

df_for_plot = merge(df_for_plot,exp_df[,c('DepMap_ID','MDM2','RB1','E2F1')])

ggscatter(df_for_plot%>%filter(TP53!='damaging'),
        x='MDM2',y='USP7',title='tp53 no damaging',
        add = "reg.line",
        conf.int = TRUE, # Add confidence interval
        cor.coef = TRUE,
        cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n"))
ggscatter(df_for_plot%>%filter(TP53=='damaging'),
        x='MDM2',y='USP7',title='tp53 damaging',
        add = "reg.line",
        conf.int = TRUE, # Add confidence interval
        cor.coef = TRUE,
        cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n"))


# 综合筛选突变对CERES的影响










