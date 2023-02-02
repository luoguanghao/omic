import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

%load_ext rpy2.ipython
%R library(ggplot2)
%R library(ggthemes)
%R library(ggbeeswarm)
%R library(dplyr)
############
# CRISPR
############
annotation_df = pd.read_csv('../data/sample_info.csv')
CRISPR_df = pd.read_csv('../data/CRISPR_gene_effect.csv')

CRISPR_df.columns = [i.split(' ')[0] for i in CRISPR_df.columns]

### all cell line ###
gene = 'CSNK1A1'

melt_df = pd.concat( [
    annotation_df.set_index('DepMap_ID').loc[set(CRISPR_df['DepMap_ID'])&set(annotation_df['DepMap_ID'])][['lineage_subtype']],
    CRISPR_df.set_index('DepMap_ID').loc[set(CRISPR_df['DepMap_ID'])&set(annotation_df['DepMap_ID'])][[gene]]
    ],axis=1 )
melt_df = melt_df.loc[pd.isna(melt_df['lineage_subtype'])==False]
melt_df.columns = ['lineage_subtype','CERES']

#melt_df = melt_df.loc[set(TP53_status_df.index)&set(melt_df.index)]
#melt_df['TP53_status'] = TP53_status_df.loc[melt_df.index]

##



## plot
# plot_melt_df = plot_melt_df.loc[plot_melt_df['TP53_status']=='WT']
plot_melt_df = melt_df

count_df = pd.concat([plot_melt_df.groupby('lineage_subtype').mean(),
                    plot_melt_df.groupby('lineage_subtype').count()],axis=1)
count_df = count_df.iloc[:,0:2]
count_df.columns = ['mean','sample_count']
ct = np.array(count_df.loc[count_df['sample_count']>0].index)

plot_melt_df = plot_melt_df.loc[ [i in ct for i in plot_melt_df['lineage_subtype']] ]

%R -i plot_melt_df

%R melt_df_mean=group_by(plot_melt_df, lineage_subtype) %>% summarize_each(funs(mean))
%R plot_melt_df$lineage_subtype = factor(plot_melt_df$lineage_subtype, levels = as.vector(melt_df_mean[order(melt_df_mean$CERES),][['lineage_subtype']]))

%R p = ggplot(plot_melt_df,aes(x = lineage_subtype, y = CERES,fill=lineage_subtype)) + geom_violin(color='white')
%R p = p + theme(text=element_text(size=25,  family="Times New Roman"),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
%R p = p + geom_jitter(aes(),fill='white',position = position_jitter(width = 0.15),shape=21)
%R p = p + stat_summary(fun.y=mean, geom="point")
%R p = p + guides(fill=FALSE)
%R p = p + ylab('CRISPR CERES score') + labs(title='CSNK1A1 CRISPR CERES in all cell lines')
%R p = p + geom_hline(aes(yintercept=mean(plot_melt_df$CERES)),linetype="dashed")
#%R plot(p)
%R ggsave("CSNK1A1_CERES_1.png", units="in", dpi=70, width=20, height=10, device="png")


### only tp53 ###
## add mut info
mut_df = pd.read_csv('../data/CCLE_mutations.csv')
mut_df.loc[mut_df['Hugo_Symbol']=='TP53'].to_csv('TP53_ccle.tsv',sep='\t')
TP53_WT_sample = list(set(mut_df['DepMap_ID'])-set(mut_df.set_index('Hugo_Symbol',drop=False).loc['TP53']['DepMap_ID']))
TP53_mut_sample = list(set(mut_df.set_index('Hugo_Symbol',drop=False).loc['TP53']['DepMap_ID']))

TP53_status_df = pd.DataFrame({
    'DepMap_ID' : TP53_WT_sample + TP53_mut_sample,
    'TP53_status' :  ['WT']*len(TP53_WT_sample)+['mut']*len(TP53_mut_sample)
}).set_index('DepMap_ID')

###
gene = 'CSNK1A1'

melt_df = pd.concat( [
    annotation_df.set_index('DepMap_ID').loc[set(CRISPR_df['DepMap_ID'])&set(annotation_df['DepMap_ID'])][['lineage_subtype']],
    CRISPR_df.set_index('DepMap_ID').loc[set(CRISPR_df['DepMap_ID'])&set(annotation_df['DepMap_ID'])][[gene]]
    ],axis=1 )
melt_df = melt_df.loc[pd.isna(melt_df['lineage_subtype'])==False]
melt_df.columns = ['lineage_subtype','CERES']

melt_df = melt_df.loc[set(TP53_status_df.index)&set(melt_df.index)]
melt_df['TP53_status'] = TP53_status_df.loc[melt_df.index]

##

## plot WT
plot_melt_df = melt_df.loc[melt_df['TP53_status']=='WT']

count_df = pd.concat([plot_melt_df.groupby('lineage_subtype').mean(),
                    plot_melt_df.groupby('lineage_subtype').count()],axis=1)
count_df = count_df.iloc[:,0:2]
count_df.columns = ['mean','sample_count']
ct = np.array(count_df.loc[count_df['sample_count']>0].index)

plot_melt_df = plot_melt_df.loc[ [i in ct for i in plot_melt_df['lineage_subtype']] ]

%R -i plot_melt_df

%R melt_df_mean=group_by(plot_melt_df, lineage_subtype) %>% summarize_each(funs(mean))
%R plot_melt_df$lineage_subtype = factor(plot_melt_df$lineage_subtype, levels = as.vector(melt_df_mean[order(melt_df_mean$CERES),][['lineage_subtype']]))

%R p = ggplot(plot_melt_df,aes(x = lineage_subtype, y = CERES,fill=lineage_subtype)) + geom_violin(color='white')
%R p = p + theme(text=element_text(size=25,  family="Times New Roman"),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
%R p = p + geom_jitter(aes(),fill='white',position = position_jitter(width = 0.15),shape=21)
%R p = p + stat_summary(fun.y=mean, geom="point")
%R p = p + guides(fill=FALSE)
%R p = p + ylab('CRISPR CERES score') + labs(title='CSNK1A1 CRISPR CERES in TP53_wt cell lines')
%R p = p + geom_hline(aes(yintercept=mean(plot_melt_df$CERES)),linetype="dashed")
#%R plot(p)
%R ggsave("CSNK1A1_CERES_tp53WT_1.png", units="in", dpi=70, width=20, height=10, device="png")



## plot mut
'''
gene = 'CSNK1A1'

melt_df = pd.concat( [
    annotation_df.set_index('DepMap_ID').loc[set(CRISPR_df['DepMap_ID'])&set(annotation_df['DepMap_ID'])][['lineage_subtype']],
    CRISPR_df.set_index('DepMap_ID').loc[set(CRISPR_df['DepMap_ID'])&set(annotation_df['DepMap_ID'])][[gene]]
    ],axis=1 )
melt_df = melt_df.loc[pd.isna(melt_df['lineage_subtype'])==False]
melt_df.columns = ['lineage_subtype','CERES']

melt_df = melt_df.loc[set(TP53_status_df.index)&set(melt_df.index)]
melt_df['TP53_status'] = TP53_status_df.loc[melt_df.index]
'''


plot_melt_df = melt_df.loc[melt_df['TP53_status']=='mut']

count_df = pd.concat([plot_melt_df.groupby('lineage_subtype').mean(),
                    plot_melt_df.groupby('lineage_subtype').count()],axis=1)
count_df = count_df.iloc[:,0:2]
count_df.columns = ['mean','sample_count']
ct = np.array(count_df.loc[count_df['sample_count']>0].index)

plot_melt_df = plot_melt_df.loc[ [i in ct for i in plot_melt_df['lineage_subtype']] ]

%R -i plot_melt_df

%R melt_df_mean=group_by(plot_melt_df, lineage_subtype) %>% summarize_each(funs(mean))
%R plot_melt_df$lineage_subtype = factor(plot_melt_df$lineage_subtype, levels = as.vector(melt_df_mean[order(melt_df_mean$CERES),][['lineage_subtype']]))

%R p = ggplot(plot_melt_df,aes(x = lineage_subtype, y = CERES,fill=lineage_subtype)) + geom_violin(color='white')
%R p = p + theme(text=element_text(size=25,  family="Times New Roman"),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
%R p = p + geom_jitter(fill='white',position = position_jitter(width = 0.15),shape=21)
%R p = p + stat_summary(fun.y=mean, geom="point")
%R p = p + guides(fill=FALSE)
%R p = p + ylab('CRISPR CERES score') + labs(title='CSNK1A1 CRISPR CERES in TP53_mut cell lines')
%R p = p + geom_hline(aes(yintercept=mean(plot_melt_df$CERES)),linetype="dashed")
#%R plot(p)
%R ggsave("CSNK1A1_CERES_tp53WT_2.png", units="in", dpi=70, width=20, height=10, device="png")













############
# RNAi
############
















