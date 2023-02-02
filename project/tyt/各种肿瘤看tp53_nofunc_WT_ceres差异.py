




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


mut_df_anno_profun_ctype_df = pd.read_csv('~/my_project/tp53_function/mut_df_anno_profun_ctype_df.tsv', sep='\t',index_col=0)


## data #
annotation_df = pd.read_csv('../data/sample_info.csv')
mut_df = pd.read_csv('../data/CCLE_mutations.csv')
CRISPR_df = pd.read_csv('../data/CRISPR_gene_effect.csv')
##

# clean data #
CRISPR_df.columns = [i.split(' ')[0] for i in CRISPR_df.columns]
# pd.merge(annotation_df,CRISPR_df,on='DepMap_ID',how='inner') # 也许应该先给crispr矩阵注释上cell line信息
#anno_CRISPR_df = pd.merge(
#    annotation_df[['DepMap_ID','CCLE_Name','lineage','lineage_subtype','lineage_sub_subtype']],
#    CRISPR_df,on='DepMap_ID',how='inner')
#


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

## Add TA information: there is a rule
cl_TA_ls = []
for cl in list(melt_df.index):
# for cl in list(test_cl):
    if cl in TP53_WT_sample:
        cl_TA = 'WT'
    else:
        cl_TA = '0'
        try:
            TAs = list(mut_df_anno_profun_ctype_df.set_index('DepMap_ID').loc[[cl]]['TransactivationClass'])
            # print(TAs)
            if 'non-functional' in TAs:
                cl_TA = 'non-functional'
            elif 'partially functional' in TAs:
                cl_TA = 'partially functional'
            elif 'functional' in TAs:
                cl_TA = 'functional'
        except:
            cl_TA = '0'
    cl_TA_ls.append(cl_TA)
melt_df['TP53_TA'] = cl_TA_ls



#######################
## 对各个lineage进行计算
#######################
lineage_subtype_ls = set(melt_df['lineage_subtype'] )

# compare_means ######
%R -i melt_df
%R res_vec = list(lineage_subtype=c(),mean_wt=c(),mean_nfuncmut=c(),pvalue=c(),adjpvalue=c(),method=c()) # non-functional mut
for ctype in lineage_subtype_ls:
    #try:
    %R -i ctype
    %R sub_melt_df = melt_df %>% filter((lineage_subtype==ctype)&(TP53_TA%in%c('non-functional','WT') ))
    levels = %R unique(sub_melt_df$TP53_TA)
    # print(levels)
    if len(levels)>1:
        %R mean_cal_res = sub_melt_df %>% group_by(TP53_TA) %>% summarise(m = mean(CERES))
        %R comp_res = compare_means(CERES~TP53_TA, data=sub_melt_df)
        %R res_vec$lineage_subtype = c(res_vec$lineage_subtype,ctype)
        %R res_vec$mean_nfuncmut = c(res_vec$mean_nfuncmut,mean_cal_res[['m']][1]) # non-functional mut
        %R res_vec$mean_wt = c(res_vec$mean_wt,mean_cal_res[['m']][2])  
        %R res_vec$pvalue = c(res_vec$pvalue, comp_res[['p']])
        %R res_vec$adjpvalue = c(res_vec$adjpvalue, comp_res[['p.adj']])
        %R res_vec$method = c(res_vec$method, comp_res[['method']])
    elif len(levels)==1:
        if levels[0]=='WT':
            %R mean_cal_res = sub_melt_df %>% group_by(TP53_TA) %>% summarise(m = mean(CERES))

            %R res_vec$lineage_subtype = c(res_vec$lineage_subtype,ctype)
            %R res_vec$mean_nfuncmut = c(res_vec$mean_nfuncmut,'888')
            %R res_vec$mean_wt = c(res_vec$mean_wt,mean_cal_res[['m']][1])  
            %R res_vec$pvalue = c(res_vec$pvalue, '888')
            %R res_vec$adjpvalue = c(res_vec$adjpvalue, '888')
            %R res_vec$method = c(res_vec$method, '-')           
        elif levels[0]=='mut':
            %R mean_cal_res = sub_melt_df %>% group_by(TP53_TA) %>% summarise(m = mean(CERES))

            %R res_vec$lineage_subtype = c(res_vec$lineage_subtype,ctype)
            %R res_vec$mean_nfuncmut = c(res_vec$mean_nfuncmut,mean_cal_res[['m']][1])
            %R res_vec$mean_wt = c(res_vec$mean_wt,'888')  
            %R res_vec$pvalue = c(res_vec$pvalue, '888')
            %R res_vec$adjpvalue = c(res_vec$adjpvalue, '888')
            %R res_vec$method = c(res_vec$method, '-')              
    else:
        continue
res_df = %R as.data.frame(res_vec)
res_df.sort_values('pvalue')
















