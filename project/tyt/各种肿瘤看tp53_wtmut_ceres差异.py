'''
新的对cell line的注释方法，有注释的就注释，没有的就标注na
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

%load_ext rpy2.ipython
%R library(ggplot2)
%R library(ggthemes)
%R library(ggbeeswarm)
%R library(dplyr)
%R library(ggpubr)

## data #
annotation_df = pd.read_csv('../data/sample_info.csv')
mut_df = pd.read_csv('../data/CCLE_mutations.csv')
CRISPR_df = pd.read_csv('../data/CRISPR_gene_effect.csv')
##

# clean data #
CRISPR_df.columns = [i.split(' ')[0] for i in CRISPR_df.columns]

# get wt/mut sample #
# mut_df.loc[mut_df['Hugo_Symbol']=='TP53'].to_csv('TP53_ccle.tsv',sep='\t')
TP53_WT_sample = list(set(mut_df['DepMap_ID'])-set(mut_df.set_index('Hugo_Symbol',drop=False).loc['TP53']['DepMap_ID']))
TP53_mut_sample = list(set(mut_df.set_index('Hugo_Symbol',drop=False).loc['TP53']['DepMap_ID']))

TP53_status_df = pd.DataFrame({
    'DepMap_ID' : TP53_WT_sample + TP53_mut_sample,
    'TP53_status' :  ['WT']*len(TP53_WT_sample)+['mut']*len(TP53_mut_sample)
}).set_index('DepMap_ID')


# sp gene #

gene = 'CSNK1A1'

# create melf df #
melt_df = pd.concat( [
    annotation_df.set_index('DepMap_ID').loc[set(CRISPR_df['DepMap_ID'])&set(annotation_df['DepMap_ID'])][['lineage_subtype']],
    CRISPR_df.set_index('DepMap_ID').loc[set(CRISPR_df['DepMap_ID'])&set(annotation_df['DepMap_ID'])][[gene]]
    ],axis=1 )
melt_df = melt_df.loc[pd.isna(melt_df['lineage_subtype'])==False] # 丢掉无亚型注释的
melt_df.columns = ['lineage_subtype','CERES']

melt_df = melt_df.loc[set(TP53_status_df.index)&set(melt_df.index)]
melt_df['TP53_status'] = TP53_status_df.loc[melt_df.index]

##

#######################
## 对各个lineage进行计算
#######################
lineage_subtype_ls = set(melt_df['lineage_subtype'] )

# t.test 舍去 ######

%R res_vec = list(lineage_subtype=c(),mean_wt=c(),mean_mut=c(),pvalue=c())

for ctype in lineage_subtype_ls:
    try:
        %R -i ctype
        %R sub_melt_df = melt_df%>%filter(lineage_subtype==ctype)
        %R res = t.test(CERES~TP53_status, data=sub_melt_df)

        %R res_vec$lineage_subtype = c(res_vec$lineage_subtype,ctype)
        %R res_vec$mean_mut = c(res_vec$mean_mut,res$estimate[1])
        %R res_vec$mean_wt = c(res_vec$mean_wt,res$estimate[2])  
        %R res_vec$pvalue = c(res_vec$pvalue, res$p.value)
    except:
        %R res_vec$lineage_subtype = c(res_vec$lineage_subtype,ctype)
        %R res_vec$mean_mut = c(res_vec$mean_mut,res$estimate[1])
        %R res_vec$mean_wt = c(res_vec$mean_wt,res$estimate[2])  
        %R res_vec$pvalue = c(res_vec$pvalue, 888)

res_df = %R as.data.frame(res_vec)
res_df.sort_values('pvalue')


# compare_means ######
%R res_vec = list(lineage_subtype=c(),mean_wt=c(),mean_mut=c(),pvalue=c(),adjpvalue=c(),method=c())
for ctype in lineage_subtype_ls:
    #try:
    %R -i ctype
    %R sub_melt_df = melt_df %>% filter(lineage_subtype==ctype)
    levels = %R unique(sub_melt_df$TP53_status)
    if len(levels)>1:
        %R mean_cal_res = sub_melt_df %>% group_by(TP53_status) %>% summarise(m = mean(CERES))
        %R comp_res = compare_means(CERES~TP53_status, data=sub_melt_df)
        %R res_vec$lineage_subtype = c(res_vec$lineage_subtype,ctype)
        %R res_vec$mean_mut = c(res_vec$mean_mut,mean_cal_res[['m']][1])
        %R res_vec$mean_wt = c(res_vec$mean_wt,mean_cal_res[['m']][2])  
        %R res_vec$pvalue = c(res_vec$pvalue, comp_res[['p']])
        %R res_vec$adjpvalue = c(res_vec$adjpvalue, comp_res[['p.adj']])
        %R res_vec$method = c(res_vec$method, comp_res[['method']])
    else:
        if levels[0]=='WT':
            %R mean_cal_res = sub_melt_df %>% group_by(TP53_status) %>% summarise(m = mean(CERES))

            %R res_vec$lineage_subtype = c(res_vec$lineage_subtype,ctype)
            %R res_vec$mean_mut = c(res_vec$mean_mut,'888')
            %R res_vec$mean_wt = c(res_vec$mean_wt,mean_cal_res[['m']][1])  
            %R res_vec$pvalue = c(res_vec$pvalue, '888')
            %R res_vec$adjpvalue = c(res_vec$adjpvalue, '888')
            %R res_vec$method = c(res_vec$method, '-')           
        elif levels[0]=='mut':
            %R mean_cal_res = sub_melt_df %>% group_by(TP53_status) %>% summarise(m = mean(CERES))

            %R res_vec$lineage_subtype = c(res_vec$lineage_subtype,ctype)
            %R res_vec$mean_mut = c(res_vec$mean_mut,mean_cal_res[['m']][1])
            %R res_vec$mean_wt = c(res_vec$mean_wt,'888')  
            %R res_vec$pvalue = c(res_vec$pvalue, '888')
            %R res_vec$adjpvalue = c(res_vec$adjpvalue, '888')
            %R res_vec$method = c(res_vec$method, '-')              

res_df = %R as.data.frame(res_vec)
res_df.sort_values('pvalue')



