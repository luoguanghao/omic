#
# 比较 c('GSPT1','CSNK1A1','IKZF1','ZFP91','IKZF3') 基因在敲除和敲低下依赖性的不同
#

library(tidyverse)
library(ggunchained)
library(ggpubr)
source('/home/lgh/my_project/omic/modules/small_tools/small_tools.R')
source('/home/lgh/my_project/omic/modules/prepare_data/depmap_data.R')



info_df = read_tsv('/home/lgh/my_project/data/depmap/cleaned/sample_info.tsv')
rnai_df = read_tsv('/home/lgh/my_project/data/depmap/cleaned/RNAi_gene_effect_cleaned.tsv')
mut_df = read_tsv('/home/lgh/my_project/data/depmap/cleaned/CCLE_mutations_cleaned.tsv')

df_for_sel_sample = data.frame(t(rnai_df%>%column_to_rownames('Hugo_Symbol')))%>%rownames_to_column('sample')
data3 <- reshape2::melt(df_for_sel_sample[,c('sample','CSNK1A1','IKZF1','ZFP91')], id.vars = c("sample"))

#######
# plot
#######

subtype = 'ALL'

celline = get_depmap_subtype_sample(info_df,mut_df=NULL,exp_df=NULL,subtype_term='lineage_subtype',subtype_name=subtype)
tp53_stat_df = get_mutation_wt_sample(mut_df,'TP53')
df_for_plot = merge(data3,tp53_stat_df)%>%filter(sample%in%celline)

p=ggviolin(df_for_plot,x='variable',y='value',fill='mutORwt',add='boxplot')+stat_compare_means(aes(group =mutORwt))
p = p+labs(title=subtype)
p = p+ylab('RNAi CERES')
p = p + geom_hline(aes(yintercept=-0.5),linetype ="dotted")
p

# 旧风格图
# df_for_volin$gene = factor(df_for_volin$gene,levels=c('GSPT1','CSNK1A1','IKZF1','ZFP91','IKZF3'))
p = ggplot(df_for_plot,aes(x = variable, y = value,fill=variable)) + geom_violin(color='white')
p=p+geom_jitter(fill='white',width =0.1,shape = 21,size=2.5)
# p=p + geom_boxplot(color='white',size=.5, width= .1)
p = p+stat_summary(fun.y=mean, geom="point")
p = p + geom_hline(aes(yintercept=-0.5),linetype ="dotted")
p = p+labs(title="ALL")
p = p+ylab('RNAi CERES')
p=p+theme(text=element_text(size=16,  family="Arial", face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5),panel.background=element_blank(),,axis.line = element_line(colour="black"))
# ,panel.background=element_blank(),,axis.line = element_line(colour="black"),panel.grid=element_line(colour="grey")
print(p)

### ========
subtype = 'AML'

celline = get_depmap_subtype_sample(info_df,mut_df=NULL,exp_df=NULL,subtype_term='lineage_subtype',subtype_name=subtype)
tp53_stat_df = get_mutation_wt_sample(mut_df,'TP53')
df_for_plot = merge(data3,tp53_stat_df)%>%filter(sample%in%celline)

p=ggviolin(df_for_plot,x='variable',y='value',fill='mutORwt',add='boxplot')+stat_compare_means(aes(group =mutORwt))
p = p+labs(title=subtype)
p = p+ylab('RNAi CERES')
p = p + geom_hline(aes(yintercept=-0.5),linetype ="dotted")
p


# 旧风格图
# df_for_volin$gene = factor(df_for_volin$gene,levels=c('GSPT1','CSNK1A1','IKZF1','ZFP91','IKZF3'))
p = ggplot(df_for_plot,aes(x = variable, y = value,fill=variable)) + geom_violin(color='white')
p=p+geom_jitter(fill='white',width =0.1,shape = 21,size=2.5)
# p=p + geom_boxplot(color='white',size=.5, width= .1)
p = p+stat_summary(fun.y=mean, geom="point")
p = p + geom_hline(aes(yintercept=-0.5),linetype ="dotted")
p = p+labs(title="ALL")
p = p+ylab('RNAi CERES')
p=p+theme(text=element_text(size=16,  family="Arial", face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5),panel.background=element_blank(),,axis.line = element_line(colour="black"))
# ,panel.background=element_blank(),,axis.line = element_line(colour="black"),panel.grid=element_line(colour="grey")
print(p)






##########
# CRISPR
##########





info_df = read_tsv('/home/lgh/my_project/data/depmap/cleaned/sample_info.tsv')

crispr_df = read_tsv('/home/lgh/my_project/data/depmap/cleaned/CRISPR_gene_effect_cleaned_2.tsv')

mut_df = read_tsv('/home/lgh/my_project/data/depmap/cleaned/CCLE_mutations_cleaned.tsv')



df_for_sel_sample = data.frame(t(crispr_df%>%column_to_rownames('Hugo_Symbol')))%>%rownames_to_column('sample')

data3 <- reshape2::melt(df_for_sel_sample[,c('sample','CSNK1A1','IKZF1','ZFP91')], id.vars = c("sample"))

source('/home/lgh/my_project/omic/modules/small_tools/small_tools.R')
source('/home/lgh/my_project/omic/modules/prepare_data/depmap_data.R')

subtype = 'ALL'

celline = get_depmap_subtype_sample(info_df,mut_df=NULL,exp_df=NULL,subtype_term='lineage_subtype',subtype_name=subtype)
tp53_stat_df = get_mutation_wt_sample(mut_df,'TP53')
df_for_plot = merge(data3,tp53_stat_df)%>%filter(sample%in%celline)

p=ggviolin(df_for_plot,x='variable',y='value',fill='mutORwt',add='boxplot')+stat_compare_means(aes(group =mutORwt))
p = p+labs(title=subtype)
p = p+ylab('CRISPR CERES')
p = p + geom_hline(aes(yintercept=-0.5),linetype ="dotted")
p

# 旧风格图
# df_for_volin$gene = factor(df_for_volin$gene,levels=c('GSPT1','CSNK1A1','IKZF1','ZFP91','IKZF3'))
p = ggplot(df_for_plot,aes(x = variable, y = value,fill=variable)) + geom_violin(color='white')
p=p+geom_jitter(fill='white',width =0.1,shape = 21,size=2.5)
# p=p + geom_boxplot(color='white',size=.5, width= .1)
p = p+stat_summary(fun.y=mean, geom="point")
p = p + geom_hline(aes(yintercept=-0.5),linetype ="dotted")
p = p+labs(title="ALL")
p = p+ylab('RNAi CERES')
p=p+theme(text=element_text(size=16,  family="Arial", face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5),panel.background=element_blank(),,axis.line = element_line(colour="black"))
# ,panel.background=element_blank(),,axis.line = element_line(colour="black"),panel.grid=element_line(colour="grey")
print(p)

## 
##


subtype = 'AML'

celline = get_depmap_subtype_sample(info_df,mut_df=NULL,exp_df=NULL,subtype_term='lineage_subtype',subtype_name=subtype)
tp53_stat_df = get_mutation_wt_sample(mut_df,'TP53')
df_for_plot = merge(data3,tp53_stat_df)%>%filter(sample%in%celline)

p=ggviolin(df_for_plot,x='variable',y='value',fill='mutORwt',add='boxplot')+stat_compare_means(aes(group =mutORwt))
p = p+labs(title=subtype)
p = p+ylab('CRISPR CERES')
p = p + geom_hline(aes(yintercept=-0.5),linetype ="dotted")
p

# 旧风格图
# df_for_volin$gene = factor(df_for_volin$gene,levels=c('GSPT1','CSNK1A1','IKZF1','ZFP91','IKZF3'))
p = ggplot(df_for_plot,aes(x = variable, y = value,fill=variable)) + geom_violin(color='white')
p=p+geom_jitter(fill='white',width =0.1,shape = 21,size=2.5)
# p=p + geom_boxplot(color='white',size=.5, width= .1)
p = p+stat_summary(fun.y=mean, geom="point")
p = p + geom_hline(aes(yintercept=-0.5),linetype ="dotted")
p = p+labs(title="ALL")
p = p+ylab('RNAi CERES')
p=p+theme(text=element_text(size=16,  family="Arial", face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5),panel.background=element_blank(),,axis.line = element_line(colour="black"))
# ,panel.background=element_blank(),,axis.line = element_line(colour="black"),panel.grid=element_line(colour="grey")
print(p)





















