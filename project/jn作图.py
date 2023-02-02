# http://n102.yfish.x:8888/notebooks/top/y/Archive/Bondi/jupyter/jianan/5-10%E5%88%86%E6%9E%90%20%E4%BB%BF%E7%85%A7AML%E7%9A%84%E5%88%86%E6%9E%90/plot2%20boxplot.ipynb

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

from rpy2.robjects.conversion import localconverter

from sklearn import metrics

import itertools
import warnings

import math

import pandas as pd
import pylab
import numpy as np


from sklearn.linear_model import enet_path
from sklearn import preprocessing
from sklearn import model_selection
from sklearn import linear_model # must use the module rather than classes to

import matplotlib.pyplot as plt

import random
%load_ext rpy2.ipython
%R library(ggplot2)
%R library(ggthemes)
#%R library(ggbeeswarm)
%R library(dplyr)

import pandas as pd
import numpy as np
import math, csv, pickle
import random, os, re
from sklearn.linear_model import ElasticNetCV
from sklearn.datasets import make_regression

from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression




%R library(RColorBrewer)

%R mycolors = c('#FF0000','#225EA8')



def ggplot2_boxplot(df, x_level, y_value, stat_method=np.nan, ref_group = ".all.", title=None, pairwise=False):
    '''
    input:
        df: 2 columns: value, level
            ! 注意，level列应该为factor

        stat_method: "t.test" or np.nan or other...
    '''
    if title is None:
        filename = 'plot.png'
        title = 'Plot: %s ~ %s'%(y_value,x_level)
    else:
        filename = '%s.png'%title
        title = '%s: %s ~ %s'%(title,y_value,x_level)
    %R library(ggplot2)
    %R library(ggpubr)        
    %R -i df
    %R -i x_level
    %R -i y_value
    %R -i stat_method
    %R if (is.nan(stat_method)){stat_method=NULL}
    %R -i ref_group
    %R -i title
    # %R p <- ggboxplot(df, x = x_level, y = y_value,color = x_level, add = "jitter")
    %R p <- ggboxplot(df, x = x_level, y = y_value,color = x_level)
    %R p = p+labs(title=title)
    #if stat_method is not np.nan:
    if pairwise:
        y_hline = np.mean(df[y_value])
        %R -i y_hline
        %R p = p + stat_compare_means(method=stat_method,label = "p.signif",ref.group=ref_group,hide.ns = FALSE)+geom_hline(yintercept = y_hline, linetype = 2)
        # %R plot(p + stat_compare_means(method=stat_method,label = "p.signif",ref.group=ref_group,hide.ns = FALSE)
        #%R plot(p+geom_hline(yintercept = y_hline, linetype = 2))
    else:
        %R p = p + stat_compare_means(method=stat_method)
        %R fo = as.formula(paste(y_value,x_level,sep='~'))
        r_stat = %R stat = compare_means(fo,data=df,method=stat_method)
        p_adj = %R stat$p.adj
        method = %R stat$method
    
    %R print(p)
    #%R -i filename
    #%R suppressMessages(ggsave(filename,p))
    #%R png(filename,width = 5,height = 10)
    #%R plot(p)
    #%R dev.off()
    #else:
    #    %R plot(p + stat_compare_means())
    'ggplot2!'
    pl = %R p
    return {'r_stat':r_stat,'p_adj':p_adj[0],'method':method[0]}


from collections import Counter
def get_common_mut(Variants_df, n=10):
    sample_count = len(set(Variants_df['Tumor_Sample_Barcode']))
    
    mut_count_df = pd.DataFrame(
        Counter(
            Variants_df['Hugo_Symbol']).items(),columns=['gene','count']).set_index('gene').sort_values('count',ascending=False)    
    
    mut_count_df['frequency'] = mut_count_df['count']/sample_count
    #most_common_res = Counter(Variants_df['Hugo_Symbol']).most_common(10)
    #top_genes = [g[0] for g in most_common_res]
    top_gene_ls = list(mut_count_df.index[:n])
    return {'mut_count_df':mut_count_df,
           'top_gene_ls':top_gene_ls}




# load vizome data
# clean data
# drug sen; mut

# import data
drug_sen_path = '/y/Archive/Bondi/data/VIZOME/nature_aml_drug_sen.tsv'
drug_sen_df = pd.read_csv(drug_sen_path, sep='\t', index_col=0)

rpkm_path = '/y/Archive/Bondi/data/VIZOME/nature_aml_log2_fpkm.txt'
rpkm_expr_mat = pd.read_csv(rpkm_path, sep='\t', index_col=0)

Variants_path =  '/y/Archive/Bondi/data/VIZOME/LAML_nature_PDC.maf'
Variants_df = pd.read_csv(Variants_path, sep='\t')

all_drug = list(set(drug_sen_df.index))

################################3

FLT3_unsen_drug_ls = ['Panobinostat']

CEBPA_unsen_drug_ls = ['GDC-0941', 'BEZ235', 'PI-103', 
           'INK-128', '17-AAG (Tanespimycin)','Go6976']

TP53_unsen_d_l = ['Lenvatinib','Regorafenib (BAY 73-4506)', 
                  'Cabozantinib','Sorafenib', 'Tivozanib (AV-951)','Pelitinib (EKB-569)']

NRAS_KRAS_unsen_d_l = ['Quizartinib (AC220)', 'Dovitinib (CHIR-258)',
                       'Lenvatinib', 'Cabozantinib', 'Vargetef', 
                       'Foretinib (XL880)', 'Tivozanib (AV-951)']

IDH_unsen_d_l = ['Bosutinib (SKI-606)', 'CYT387']

PTPN11_unsen_d_l = ['ABT-737', 'Venetoclax']

################################
# 合成res_df

res_df = None


sp_mut = 'FLT3'     
pandel = sp_mut
for drug in FLT3_unsen_drug_ls:
# for drug in set(drug_sen_df.index):
    print('.',end='')
    useful_sample_ls = list(set(Variants_df['Tumor_Sample_Barcode'])&set(drug_sen_df.loc[drug]['lab_id']))
    useful_Variants_df = Variants_df.loc[[i in useful_sample_ls for i in Variants_df['Tumor_Sample_Barcode']]]

    # for sp_mut in get_common_mut(useful_Variants_df)['mut_count_df'].iloc[:].index:
    sp_Variants_df = useful_Variants_df.loc[useful_Variants_df['Hugo_Symbol']==sp_mut]
    sp_sample_ls = list(set(sp_Variants_df['Tumor_Sample_Barcode']))
    unsp_sample_ls = list(set(useful_sample_ls)-set(sp_sample_ls))

    if len(sp_sample_ls)<3:
        print('Drug %s finish, the last mutatioin is %s, the sample_num is %s'%(drug,sp_mut,len(sp_sample_ls)))
        break

    sp_mut_drug_sen_df = drug_sen_df.loc[drug].set_index('lab_id').loc[sp_sample_ls][['ic50']]
    unsp_mut_drug_sen_df = drug_sen_df.loc[drug].set_index('lab_id').loc[unsp_sample_ls][['ic50']]

    unsp_mut_drug_sen_df['flag']='WT'
    sp_mut_drug_sen_df['flag']='Mutation'

    df_for_cal = pd.concat([sp_mut_drug_sen_df,unsp_mut_drug_sen_df],axis=0)
    # df_for_cal['lab_id'] = df_for_cal.index
    df_for_cal['drug'] = drug
    df_for_cal['pandel'] = pandel

    if res_df is None:
        res_df = df_for_cal
    else:
        res_df = pd.concat([res_df,df_for_cal],axis=0)

#############################################

sp_mut = 'CEBPA'     
pandel = sp_mut
for drug in CEBPA_unsen_drug_ls:
# for drug in set(drug_sen_df.index):
    print('.',end='')
    useful_sample_ls = list(set(Variants_df['Tumor_Sample_Barcode'])&set(drug_sen_df.loc[drug]['lab_id']))
    useful_Variants_df = Variants_df.loc[[i in useful_sample_ls for i in Variants_df['Tumor_Sample_Barcode']]]

    # for sp_mut in get_common_mut(useful_Variants_df)['mut_count_df'].iloc[:].index:
    sp_Variants_df = useful_Variants_df.loc[useful_Variants_df['Hugo_Symbol']==sp_mut]
    sp_sample_ls = list(set(sp_Variants_df['Tumor_Sample_Barcode']))
    unsp_sample_ls = list(set(useful_sample_ls)-set(sp_sample_ls))

    if len(sp_sample_ls)<3:
        print('Drug %s finish, the last mutatioin is %s, the sample_num is %s'%(drug,sp_mut,len(sp_sample_ls)))
        break

    sp_mut_drug_sen_df = drug_sen_df.loc[drug].set_index('lab_id').loc[sp_sample_ls][['ic50']]
    unsp_mut_drug_sen_df = drug_sen_df.loc[drug].set_index('lab_id').loc[unsp_sample_ls][['ic50']]

    unsp_mut_drug_sen_df['flag']='WT'
    sp_mut_drug_sen_df['flag']='Mutation'

    df_for_cal = pd.concat([sp_mut_drug_sen_df,unsp_mut_drug_sen_df],axis=0)
    # df_for_cal['lab_id'] = df_for_cal.index
    df_for_cal['drug'] = drug
    df_for_cal['pandel'] = pandel

    if res_df is None:
        res_df = df_for_cal
    else:
        res_df = pd.concat([res_df,df_for_cal],axis=0)

#############################################
# res_df = pd.DataFrame(res_dict)
sp_mut = 'TP53'
pandel = 'TP53'
for drug in TP53_unsen_d_l:
# for drug in set(drug_sen_df.index):
    print('.',end='')
    useful_sample_ls = list(set(Variants_df['Tumor_Sample_Barcode'])&set(drug_sen_df.loc[drug]['lab_id']))
    useful_Variants_df = Variants_df.loc[[i in useful_sample_ls for i in Variants_df['Tumor_Sample_Barcode']]]

    # for sp_mut in get_common_mut(useful_Variants_df)['mut_count_df'].iloc[:].index:
    sp_Variants_df = useful_Variants_df.loc[useful_Variants_df['Hugo_Symbol']==sp_mut]
    sp_sample_ls = list(set(sp_Variants_df['Tumor_Sample_Barcode']))
    unsp_sample_ls = list(set(useful_sample_ls)-set(sp_sample_ls))

    if len(sp_sample_ls)<3:
        print('Drug %s finish, the last mutatioin is %s, the sample_num is %s'%(drug,sp_mut,len(sp_sample_ls)))
        break

    sp_mut_drug_sen_df = drug_sen_df.loc[drug].set_index('lab_id').loc[sp_sample_ls][['ic50']]
    unsp_mut_drug_sen_df = drug_sen_df.loc[drug].set_index('lab_id').loc[unsp_sample_ls][['ic50']]

    unsp_mut_drug_sen_df['flag']='WT'
    sp_mut_drug_sen_df['flag']='Mutation'

    df_for_cal = pd.concat([sp_mut_drug_sen_df,unsp_mut_drug_sen_df],axis=0)
    # df_for_cal['lab_id'] = df_for_cal.index
    df_for_cal['drug'] = drug
    df_for_cal['pandel'] = pandel

    if res_df is None:
        res_df = df_for_cal
    else:
        res_df = pd.concat([res_df,df_for_cal],axis=0)

#############################################
sp_mut = 'KRAS'  ##### N/KRAS
pandel = 'N/KRAS'
for drug in NRAS_KRAS_unsen_d_l:
# for drug in set(drug_sen_df.index):
    print('.',end='')
    useful_sample_ls = list(set(Variants_df['Tumor_Sample_Barcode'])&set(drug_sen_df.loc[drug]['lab_id']))
    useful_Variants_df = Variants_df.loc[[i in useful_sample_ls for i in Variants_df['Tumor_Sample_Barcode']]]

    # for sp_mut in get_common_mut(useful_Variants_df)['mut_count_df'].iloc[:].index:
    sp_Variants_df = useful_Variants_df.loc[ 
        useful_Variants_df['Hugo_Symbol'].apply(lambda x: pd.isna(re.match('^(N|K)RAS',x))==False ) ]
    sp_sample_ls = list(set(sp_Variants_df['Tumor_Sample_Barcode']))
    unsp_sample_ls = list(set(useful_sample_ls)-set(sp_sample_ls))

    if len(sp_sample_ls)<3:
        print('Drug %s finish, the last mutatioin is %s, the sample_num is %s'%(drug,sp_mut,len(sp_sample_ls)))
        break

    sp_mut_drug_sen_df = drug_sen_df.loc[drug].set_index('lab_id').loc[sp_sample_ls][['ic50']]
    unsp_mut_drug_sen_df = drug_sen_df.loc[drug].set_index('lab_id').loc[unsp_sample_ls][['ic50']]

    unsp_mut_drug_sen_df['flag']='WT'
    sp_mut_drug_sen_df['flag']='Mutation'

    df_for_cal = pd.concat([sp_mut_drug_sen_df,unsp_mut_drug_sen_df],axis=0)
    # df_for_cal['lab_id'] = df_for_cal.index
    df_for_cal['drug'] = drug
    df_for_cal['pandel'] = pandel

    if res_df is None:
        res_df = df_for_cal
    else:
        res_df = pd.concat([res_df,df_for_cal],axis=0)

#############################################
sp_mut = 'IDH'
pandel = 'IDH'
for drug in IDH_unsen_d_l:
# for drug in set(drug_sen_df.index):
    print('.',end='')
    useful_sample_ls = list(set(Variants_df['Tumor_Sample_Barcode'])&set(drug_sen_df.loc[drug]['lab_id']))
    useful_Variants_df = Variants_df.loc[[i in useful_sample_ls for i in Variants_df['Tumor_Sample_Barcode']]]

    # for sp_mut in get_common_mut(useful_Variants_df)['mut_count_df'].iloc[:].index:
    sp_Variants_df = useful_Variants_df.loc[useful_Variants_df['Hugo_Symbol']==sp_mut]
    sp_sample_ls = list(set(sp_Variants_df['Tumor_Sample_Barcode']))
    unsp_sample_ls = list(set(useful_sample_ls)-set(sp_sample_ls))

    if len(sp_sample_ls)<3:
        print('Drug %s finish, the last mutatioin is %s, the sample_num is %s'%(drug,sp_mut,len(sp_sample_ls)))
        break

    sp_mut_drug_sen_df = drug_sen_df.loc[drug].set_index('lab_id').loc[sp_sample_ls][['ic50']]
    unsp_mut_drug_sen_df = drug_sen_df.loc[drug].set_index('lab_id').loc[unsp_sample_ls][['ic50']]

    unsp_mut_drug_sen_df['flag']='WT'
    sp_mut_drug_sen_df['flag']='Mutation'

    df_for_cal = pd.concat([sp_mut_drug_sen_df,unsp_mut_drug_sen_df],axis=0)
    # df_for_cal['lab_id'] = df_for_cal.index
    df_for_cal['drug'] = drug
    df_for_cal['pandel'] = pandel

    if res_df is None:
        res_df = df_for_cal
    else:
        res_df = pd.concat([res_df,df_for_cal],axis=0)

#############################################
sp_mut = 'PTPN11'
pandel = 'PTPN11'
for drug in PTPN11_unsen_d_l:
# for drug in set(drug_sen_df.index):
    print('.',end='')
    useful_sample_ls = list(set(Variants_df['Tumor_Sample_Barcode'])&set(drug_sen_df.loc[drug]['lab_id']))
    useful_Variants_df = Variants_df.loc[[i in useful_sample_ls for i in Variants_df['Tumor_Sample_Barcode']]]

    # for sp_mut in get_common_mut(useful_Variants_df)['mut_count_df'].iloc[:].index:
    sp_Variants_df = useful_Variants_df.loc[useful_Variants_df['Hugo_Symbol']==sp_mut]
    sp_sample_ls = list(set(sp_Variants_df['Tumor_Sample_Barcode']))
    unsp_sample_ls = list(set(useful_sample_ls)-set(sp_sample_ls))

    if len(sp_sample_ls)<3:
        print('Drug %s finish, the last mutatioin is %s, the sample_num is %s'%(drug,sp_mut,len(sp_sample_ls)))
        break

    sp_mut_drug_sen_df = drug_sen_df.loc[drug].set_index('lab_id').loc[sp_sample_ls][['ic50']]
    unsp_mut_drug_sen_df = drug_sen_df.loc[drug].set_index('lab_id').loc[unsp_sample_ls][['ic50']]

    unsp_mut_drug_sen_df['flag']='WT'
    sp_mut_drug_sen_df['flag']='Mutation'

    df_for_cal = pd.concat([sp_mut_drug_sen_df,unsp_mut_drug_sen_df],axis=0)
    # df_for_cal['lab_id'] = df_for_cal.index
    df_for_cal['drug'] = drug
    df_for_cal['pandel'] = pandel

    if res_df is None:
        res_df = df_for_cal
    else:
        res_df = pd.concat([res_df,df_for_cal],axis=0)


res_df['drug'] = res_df['drug'].apply(lambda x : x.split(' (')[0])

%R -i res_df
%R res_df[2:3] <- lapply( res_df[2:3], factor)
%R res_df$ic50 = res_df$ic50  ###############
%R p = ggplot(data = res_df, aes(x = drug, y = ic50, fill = factor(flag)))
#%R p = p + geom_violin(aes(color=flag))
%R p = p + geom_boxplot(outlier.colour = NA,linetype="dashed") ## 设置虚线
# %R p = p + stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fatten = NULL))  # 去掉须线的纯框
%R p = p + geom_boxplot(aes(ymin=..lower.., ymax=..upper..),fatten = NULL, outlier.size = 0.5)
%R p = p + stat_summary(fun.y=mean, geom="point",shape=8,position = position_dodge(0.8),color='black')
# %R p = p + stat_boxplot(geom = "errorbar",aes(x = drug,y=ic50,color = 'black'),width=0.2)
# %R p = p + stat_boxplot(geom = "errorbar",aes(ymax=..ymin.., fill = factor(flag)),width=0.2)
# %R p = p + coord_fixed(ratio=1.5)
%R p = p + facet_grid(.~pandel, scales = "free",space="free")

#%R plot(p+labs(title='PI3K/MTOR resistance in CEBPA mutation PDC sample')+theme(text=element_text(size=15,family="serif")))
# %R p=p+labs(title='PI3K/MTOR resistance in CEBPA mutation PDC sample')
# %R p=p+theme(text=element_text(size=15,family="serif"))
%R p=p+scale_fill_manual("flag",values = mycolors)+scale_color_manual("flag",values = mycolors)
%R p = p + theme(text=element_text(size=12,  family="Arial",face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5), panel.background = element_blank(),axis.line = element_line(colour="black"),legend.key = element_rect(fill = "white"))
%R ggsave("./res/plot2_boxplot_mut_unsen_part.png", units="in", dpi=100, width=12, height=3.5, device="png")











