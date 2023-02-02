# 有这几种图：
## 某基因在各个cell line中的dependence图
## 几个基因在某癌症中的dependence
## 几个基因在某癌症中的dependence的分布（所有基因中的分布）
## ! 针对CRISPR数据和RNAi数据 
## ! 要求这个作图框架可移植性强
#
#

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

%load_ext rpy2.ipython
%R library(ggplot2)
%R library(ggthemes)
%R library(ggbeeswarm)
%R library(dplyr)
%R library(ConsensusClusterPlus)

################
################
# CRISPR ceres #
################
################

#################################################
# 整体分布图 某个基因的表达情况，在各个细胞系的情况 #
#################################################

# 数据获取  ##########################
annotation_df = pd.read_csv('../../data/sample_info.csv')
CRISPR_df = pd.read_csv('../../data/CRISPR_gene_effect.csv')

CRISPR_df.columns = [i.split(' ')[0] for i in CRISPR_df.columns]

# 构建合并矩阵 （选择特定的基因） ##########################
gene = 'PTP4A3'

##### /////
output_name_violin = '%s_CancetType_violin_CRISPR.png'%gene

## 获取融合矩阵
melt_df = annotation_df.set_index('lineage_subtype',drop=False)[
            ['cell_line_name','stripped_cell_line_name','DepMap_ID','CCLE_Name','lineage_subtype']]

tmp_ls = set(melt_df['DepMap_ID'])&set(CRISPR_df['DepMap_ID']) # 选都有的样本
melt_df = pd.concat([melt_df.set_index('DepMap_ID',drop=False).loc[tmp_ls],
    CRISPR_df.set_index('DepMap_ID').loc[tmp_ls][[gene]]],axis=1) 


## 选取合适的lineage
count_df = pd.concat([melt_df.groupby('lineage_subtype').mean(),
                    melt_df.groupby('lineage_subtype').count()['DepMap_ID']],axis=1)
count_df.columns = ['mean','sample_count']
ct = np.array(count_df.loc[count_df['sample_count']>5].index) ## 这里可以改

## plot

plot_melt_df = melt_df.loc[[i in ct for i in melt_df['lineage_subtype']]][['lineage_subtype',gene]]
plot_melt_df.columns = ['lineage_subtype','CERES']
%R -i plot_melt_df
%R plot_melt_df = na.omit(plot_melt_df)
%R -i output_name_violin
%R -i gene
%R melt_df_mean=group_by(plot_melt_df, lineage_subtype) %>% summarize_each(funs(mean))
%R plot_melt_df$lineage_subtype = factor(plot_melt_df$lineage_subtype, levels = as.vector(melt_df_mean[order(melt_df_mean$CERES),][['lineage_subtype']]))

%R p = ggplot(plot_melt_df,aes(x = lineage_subtype, y = CERES,fill=lineage_subtype)) + geom_violin(color='white')
%R p = p + theme(text=element_text(size=15,  family="Times New Roman"),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
%R p = p + geom_jitter(aes(fill=lineage_subtype),position = position_jitter(width = 0.15),shape=21)
%R p = p + stat_summary(fun.y=mean, geom="point")
%R p = p + labs(title=gene) + ylab('CRISPR CERES score')
#%R plot(p)
%R ggsave(output_name_violin, units="in", dpi=300, width=25, height=10, device="png")

### barplot
# %R -i to_plot_boxplot_df
# %R to_plot_boxplot_df$Cell_line = factor(to_plot_boxplot_df$Cell_line, levels = as.vector(to_plot_boxplot_df[order(to_plot_boxplot_df[,3]),][['Cell_line']]))
%R p = ggplot(data=melt_df_mean, aes(x=lineage_subtype, y=CERES, fill=lineage_subtype)) + geom_bar(stat="identity")
%R p = p + theme(text=element_text(size=13,  family="Times New Roman"),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
%R ggsave("test1.png", units="in", dpi=300, width=20, height=10, device="png") # 这里画柱状图
#%R plot(p)


## 额外的分析 

### 获取前于ALL

### 计算各cell type mean score，和各cell type的sample数量合并
mean_count_df = pd.concat([melt_df.groupby('lineage_subtype').mean(),
                    melt_df.groupby('lineage_subtype').count()],axis=1)
mean_count_df.columns = ['mean','sample_count']
mean_count_df['lineage_subtype'] = mean_count_df.index

# mean_count_df.sort_values('mean').to_csv('lineage_subtype_score_mean_df.csv')

### 获取在ALL之前的cancer type的sample的score
mean_count_df = mean_count_df.sort_values('mean')
mean_count_df.index = range(mean_count_df.shape[0])

res_df = annotation_df.set_index('lineage_subtype',drop=False).loc[mean_count_df.iloc[:26]['lineage_subtype']][
            ['cell_line_name','stripped_cell_line_name','DepMap_ID','lineage_subtype']] ## 获取ALL之前的cancer type

tmp_ls = set(res_df['DepMap_ID'])&set(CRISPR_df['DepMap_ID'])
res_df = pd.concat([res_df.set_index('DepMap_ID',drop=False).loc[tmp_ls],
    CRISPR_df.set_index('DepMap_ID').loc[tmp_ls]],axis=1) ## 获取ALL之前的cancer type的cell line 把annotation_df 与 score mat合并，方便获取信息

res_df = res_df.set_index('lineage_subtype',drop=False).loc[mean_count_df.iloc[:26]['lineage_subtype']][
                    ['cell_line_name','stripped_cell_line_name','DepMap_ID','lineage_subtype',gene]]

# res_df.set_index('lineage_subtype',drop=False).loc[mean_count_df.iloc[:26]['lineage_subtype']][
#                    ['cell_line_name','stripped_cell_line_name','DepMap_ID','lineage_subtype','CSNK1A1']
#                    ].to_csv('ALL之前的肿瘤细胞株CrisprScoreList.csv')


old_cell_line_ls = open('../old_cell_line_瘤谱汇总.txt').read().strip().split('\n\n\n')

print(set(res_df['stripped_cell_line_name'])&set(old_cell_line_ls))

print(set(res_df['cell_line_name'])&set(old_cell_line_ls))

old_cell_line_ls = open('../old_cell_line_项目化合物.txt').read().strip().split('\n\n\n')

print(set(res_df['stripped_cell_line_name'])&set(old_cell_line_ls))

print(set(res_df['cell_line_name'])&set(old_cell_line_ls))






########################################
# 看目的基因在所有基因中的位置，局限于ALL #
########################################

def mylog(x):
    if x > 0.5:
        x = np.log10(x)+0.6
    if x < -0.5:
        x = -(np.log10(-x)+0.6)
    return x

## 数据
annotation_df = pd.read_csv('../../data/sample_info.csv')
CRISPR_df = pd.read_csv('../../data/CRISPR_gene_effect.csv')

CRISPR_df.columns = [i.split(' ')[0] for i in CRISPR_df.columns]

## action：

DepMap_ID_ls = []
for i_r in range(annotation_df.shape[0]):
    if annotation_df.iloc[i_r]['Subtype']==annotation_df.iloc[i_r]['Subtype']:
        if 'ALL' in annotation_df.iloc[i_r]['Subtype']:
            DepMap_ID_ls.append(annotation_df.iloc[i_r]['DepMap_ID'])

DepMap_ID_ls = list(set(DepMap_ID_ls)&set(CRISPR_df['DepMap_ID']))
select_CRISPR_df = CRISPR_df.set_index('DepMap_ID').loc[DepMap_ID_ls]
select_CRISPR_df.columns = [i.split(' ')[0] for i in select_CRISPR_df.columns]

gene_sorted_df = pd.DataFrame(np.mean(select_CRISPR_df).sort_values(),columns=['mean_value'])
# gene_sorted_df['gene'] = [i.split(' ')[0] for i in gene_sorted_df.index]
gene_sorted_df['gene'] = gene_sorted_df.index
gene_sorted_df.index = np.array(range(gene_sorted_df.shape[0]))+1
gene_sorted_df = gene_sorted_df[['gene','mean_value']]

# Action #

# gene_ls = g_l.strip().split('\n')
gene_ls = ['USP7','CSNK1A1L', 'CSNK1D', 'CSNK1E', 'CSNK1G2', 'CSNK1G3', 'CSNK1G1', 'TTBK1', 'VRK1', 'TTBK2', 'VRK2', 'CDC7', 'VRK3'][:]
vl_x_ls = [mylog(float(gene_sorted_df.loc[gene_sorted_df['gene']==g]['mean_value'])) for g in gene_ls]

## 设置上下限
y_max = max(gene_sorted_df['mean_value']) + 0.5
y_min = min(gene_sorted_df['mean_value']) - 0.5
# 设置scatter的位置
x_for_scat = list(gene_sorted_df['mean_value'] )
x_for_scat = [mylog(x) for x in x_for_scat]


# =======================================================
#y_for_scat = (1.5,0.1)
#y_for_scat = [mylog(y) for y in y_for_scat]
## 画出密度图、直方图
fig,ax=plt.subplots(figsize=(15,3))
sns.distplot(x_for_scat,rug=True,                
             #kde_kws = {"shade": False,          # 进行面积填充
             #               "color": 'darkorange',  # 设置线条颜色
             #               # 'linewidth': 1.0,     # 设置线条粗细
             #               'facecolor': 'gray'},   # 设置填充颜色
                 rug_kws = {'color': 'blue',         # 设置 rug 颜色
                            'height': 0.02})       
## 画出我们的gene的位置
for i in range(len(gene_ls)):
    ax.vlines(x=vl_x_ls[i], ymin=0, ymax=1, color='red')
    plt.text(x=vl_x_ls[i],y=1,s=gene_ls[i], ha='center', va= 'bottom',fontsize=10,rotation=90)
'''
ax.vlines(x=vl_x1, ymin=0, ymax=1, color='red')
plt.text(x=vl_x1,y=1,s=gene1, ha='center', va= 'bottom',fontsize=10,rotation=45)
ax.vlines(x=vl_x2, ymin=0, ymax=1, color='red')
plt.text(x=vl_x2,y=1,s=gene2, ha='center', va= 'bottom',fontsize=10,rotation=45)
'''
## 设置ticks,label,进行画图
up_xticks = list(np.arange(1,2,0.5))
down_xticks = list(np.arange(-2.5,-0.5,0.5))
xticks = list((np.arange(-0.5,0.5+0.05,0.1)))
xticks = [round(x,2) for x in xticks]

xtickslabel = down_xticks + xticks + up_xticks
#xtickslabel = down_xticks + [-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5] + up_xticks
xticks = [mylog(x) for x in xtickslabel]
ax.set_xticks(xticks)
ax.set_xticklabels(xtickslabel, fontsize=8)
plt.show()








###############################
# 画出各个基因在ALL的violin图 ###
###############################

## 数据
annotation_df = pd.read_csv('../../data/sample_info.csv')
CRISPR_df = pd.read_csv('../../data/CRISPR_gene_effect.csv')

CRISPR_df.columns = [i.split(' ')[0] for i in CRISPR_df.columns]

## 合并出矩阵

merge_df = annotation_df.set_index('lineage_subtype',drop=False)[
            ['cell_line_name','stripped_cell_line_name','DepMap_ID','lineage_subtype']] 

tmp_ls = set(merge_df['DepMap_ID'])&set(CRISPR_df['DepMap_ID'])
merge_df = pd.concat([merge_df.set_index('DepMap_ID',drop=False).loc[tmp_ls],
    CRISPR_df.set_index('DepMap_ID').loc[tmp_ls]],axis=1) 
# 这里依赖RNAi数据获得merge_df：
# 'cell_line_name','stripped_cell_line_name','DepMap_ID','lineage_subtype'，还有各个基因的CERES
# 接下来就可以根据这个df来做各种图
# 这里的index是 DepMap_ID


## action

gene_ls = ['USP7','CSNK1A1L', 'CSNK1D', 'CSNK1E', 'CSNK1G2', 'CSNK1G3', 
            'CSNK1G1', 'TTBK1', 'VRK1', 'TTBK2', 'VRK2', 'CDC7', 'VRK3']


## ALL
DepMap_ID_ls = []
for i_r in range(annotation_df.shape[0]):
    if annotation_df.iloc[i_r]['Subtype']==annotation_df.iloc[i_r]['Subtype']:
        if 'ALL' in annotation_df.iloc[i_r]['Subtype']:  ####### disease can be change!!
            DepMap_ID_ls.append(annotation_df.iloc[i_r]['DepMap_ID'])

df_for_volin = pd.melt(merge_df,id_vars=['cell_line_name','lineage_subtype','stripped_cell_line_name','DepMap_ID'],
        value_vars=gene_ls[:],var_name='gene',value_name='CERES_score')

tmp_ls = list(set(df_for_volin['DepMap_ID'])&set(DepMap_ID_ls))
df_for_volin = df_for_volin.set_index('DepMap_ID',drop=False).loc[tmp_ls]

%R -i df_for_volin
%R df_for_volin$'gene' = as.factor(df_for_volin$'gene')

%R datamean=group_by(df_for_volin, gene) %>% summarize_each(funs(mean))
%R df_for_volin$gene = factor(df_for_volin$gene, levels = as.vector(datamean[order(datamean$CERES_score),][['gene']]))
%R p = ggplot(df_for_volin,aes(x = gene, y = CERES_score, fill=gene)) 
%R p = p + geom_violin(color='white')
#%R p=p + geom_boxplot(color='white',fill='grey',size=.4, width= .15)
%R p=p+geom_jitter(fill='white',width =0.1,shape = 21,size=1.5)
%R p = p+stat_summary(fun.y=mean, geom="point")
%R p = p+labs(title="ALL")
%R p = p+ylab('CRISPR_score')
%R p=p+theme(text=element_text(size=16,  family="Arial"),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))

#%R p = p + geom_hline(aes(yintercept=-0.263149))
#%R print(p)
%R ggsave("test3.png", units="in", dpi=100, width=10, height=8, device="png")









################
################
# RNAi ceres   #
################
################


annotation_df = pd.read_csv('../../data/sample_info.csv')

RNAi_df = pd.read_csv('../../data/D2_combined_gene_dep_scores.csv')
RNAi_df = RNAi_df.set_index(['Unnamed: 0']).T
RNAi_df.columns = [i.split(' ')[0] for i in RNAi_df.columns]

melt_df = annotation_df.set_index('lineage_subtype',drop=False)[
            ['cell_line_name','stripped_cell_line_name','DepMap_ID','CCLE_Name','lineage_subtype']]

tmp_ls = set(melt_df['CCLE_Name'])&set(RNAi_df.index) # 选都有的样本
melt_df = pd.concat([melt_df.set_index('CCLE_Name',drop=False).loc[tmp_ls],
    RNAi_df.loc[tmp_ls][[gene]]],axis=1) 

# 这里依赖RNAi数据获得merge_df：
# 'cell_line_name','stripped_cell_line_name','DepMap_ID','CCLE_Name','lineage_subtype'，还有各个基因的CERES
# 接下来就可以根据这个df来做各种图
# 这里和CRISPR不同是cell line的名字这里是CCLE_Name，而不是depmap_ID

##### /////
output_name_violin = '%s_CancetType_violin_RNAi.png'%gene

## 选取合适的lineage
count_df = pd.concat([melt_df.groupby('lineage_subtype').mean(),
                    melt_df.groupby('lineage_subtype').count()['DepMap_ID']],axis=1)
count_df.columns = ['mean','sample_count']
ct = np.array(count_df.loc[count_df['sample_count']>5].index)

## plot

plot_melt_df = melt_df.loc[[i in ct for i in melt_df['lineage_subtype']]][['lineage_subtype',gene]]
plot_melt_df.columns = ['lineage_subtype','CERES']
%R -i plot_melt_df
%R plot_melt_df = na.omit(plot_melt_df)
%R -i output_name_violin
%R -i gene
%R melt_df_mean=group_by(plot_melt_df, lineage_subtype) %>% summarize_each(funs(mean))
%R plot_melt_df$lineage_subtype = factor(plot_melt_df$lineage_subtype, levels = as.vector(melt_df_mean[order(melt_df_mean$CERES),][['lineage_subtype']]))

%R p = ggplot(plot_melt_df,aes(x = lineage_subtype, y = CERES,fill=lineage_subtype)) + geom_violin(color='white')
%R p = p + theme(text=element_text(size=15,  family="Times New Roman"),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
%R p = p + geom_jitter(aes(fill=lineage_subtype),position = position_jitter(width = 0.15),shape=21)
%R p = p + stat_summary(fun.y=mean, geom="point")
%R p = p + labs(title=gene) + ylab('RNAi CERES score')
#%R plot(p)
%R ggsave(output_name_violin, units="in", dpi=300, width=25, height=10, device="png")










