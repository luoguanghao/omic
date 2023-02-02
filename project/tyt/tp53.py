import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

%load_ext rpy2.ipython
%R library(ggplot2)
%R library(ggthemes)
%R library(ggbeeswarm)
%R library(dplyr)

# 数据获取  ##########################
annotation_df = pd.read_csv('~/my_project/data/depmap/info/sample_info.csv')
CRISPR_df = pd.read_csv('~/my_project/data/depmap/dependence/CRISPR_gene_effect.csv')

CRISPR_df.columns = [i.split(' ')[0] for i in CRISPR_df.columns]
mut_df = pd.read_csv('~/my_project/data/depmap/genome/CCLE_mutations.csv')

############
# old but ok
############

df_for_plot = CRISPR_df.set_index('DepMap_ID')[['USP7']]
df_for_plot = df_for_plot.loc[set(df_for_plot.index) & set(mut_df['DepMap_ID'])] # 以CERES和mut_df都有的样本为底版
df_for_plot = df_for_plot.loc[pd.isna(df_for_plot['USP7'])==False]
df_for_plot['TP53_statue'] = None



dep_id_ls = []
mut_var_class_ls = []

tmp_mut_df = mut_df.set_index('Hugo_Symbol',drop=False).loc['TP53'] # 选取tp53 MUT的
for dep_id in df_for_plot.index:
    try:
        mut_var_class = tmp_mut_df.set_index('DepMap_ID',drop=False).loc[[dep_id]]['Variant_Classification']
        mut_var_class = ';'.join(set(mut_var_class))
    except:
        mut_var_class = 'WT'
    # df_for_plot.loc[dep_id]['TP53_statue'] = mut_var_class
    dep_id_ls.append(dep_id)
    mut_var_class_ls.append(mut_var_class)
    print('.',end='')
df_for_plot['TP53_statue'] = mut_var_class_ls

%R -i df_for_plot
%R df_for_plot = na.omit(df_for_plot)
%R df_for_plot$TP53_statue = factor(df_for_plot$TP53_statue)

%R df_mean=group_by(df_for_plot, TP53_statue) %>% summarize_each(funs(mean))
%R df_for_plot$TP53_statue = factor(df_for_plot$TP53_statue, levels = as.vector(df_mean[order(df_mean$USP7),][['TP53_statue']]))

%R p = ggplot(df_for_plot,aes(x = TP53_statue, y = USP7,fill=TP53_statue)) + geom_violin(color='white')
%R p = p + geom_jitter(aes(fill=TP53_statue),position = position_jitter(width = 0.1),shape=21)
%R p = p + stat_summary(fun.y=mean, geom="point")
%R p = p + theme(text=element_text(size=20,  family="Times New Roman"),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
%R ggsave('tp53mut_usp7.png', units="in", dpi=100, width=20, height=15, device="png")
# %R plot(p)


#######
# anno
#######
df_for_plot = CRISPR_df.set_index('DepMap_ID')[['USP7']]
df_for_plot = df_for_plot.loc[set(df_for_plot.index) & set(mut_df['DepMap_ID'])] # 以CERES和mut_df都有的样本为底版
df_for_plot = df_for_plot.loc[pd.isna(df_for_plot['USP7'])==False]
df_for_plot['TP53_Variant_annotation'] = None

dep_id_ls = []
mut_var_class_raw_ls = []
mut_var_class_ls = []

tmp_mut_df = mut_df.set_index('Hugo_Symbol',drop=False).loc['TP53']
for dep_id in df_for_plot.index:
    try:
        mut_var_class = list(tmp_mut_df.set_index('DepMap_ID',drop=False).loc[[dep_id]]['Variant_annotation'])
        # if len(mut_var_class)>1:
            # print(dep_id,mut_var_class)
        # mut_var_class = ';'.join(set(mut_var_class))
    except:
        mut_var_class = ['WT']
    # df_for_plot.loc[dep_id]['TP53_statue'] = mut_var_class
    mut_var_class_raw = mut_var_class
    if len(set(mut_var_class))>1:
        if 'damaging' in set(mut_var_class):
            mut_var_class = 'damaging'
        elif 'other non-conserving' in set(mut_var_class):
            mut_var_class = 'other non-conserving'
        else:
            mut_var_class = 'silent'
    else:
        mut_var_class = mut_var_class[0]
    
    dep_id_ls.append(dep_id)
    mut_var_class_ls.append(mut_var_class)
    mut_var_class_raw_ls.append(mut_var_class_raw)
    print('.',end='')
    
df_for_plot['TP53_Variant_annotation'] = list(map(str,mut_var_class_ls)) # 这里把na变成了'nan'

%R -i df_for_plot
%R df_for_plot = na.omit(df_for_plot)
%R df_for_plot$TP53_Variant_annotation = factor(df_for_plot$TP53_Variant_annotation)

%R df_mean=group_by(df_for_plot, TP53_Variant_annotation) %>% summarize_each(funs(mean))
%R df_for_plot$TP53_Variant_annotation = factor(df_for_plot$TP53_Variant_annotation, levels = as.vector(df_mean[order(df_mean$USP7),][['TP53_Variant_annotation']]))

%R p = ggplot(df_for_plot,aes(x = TP53_Variant_annotation, y = USP7,fill=TP53_Variant_annotation)) + geom_violin(color='white')
%R p = p + geom_jitter(aes(fill=TP53_Variant_annotation),position = position_jitter(width = 0.1),shape=21)
%R p = p + stat_summary(fun.y=mean, geom="point")
%R p = p + theme(text=element_text(size=30,  family="Times New Roman"),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
%R ggsave('tp53mut_usp7_factor_anno.png', units="in", dpi=100,width=10,height=10,  device="png")
# %R plot(p)






TP53_anno



















