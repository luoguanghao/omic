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

import pandas as pd
import numpy as np
import math, csv, pickle
import random, os
from sklearn.linear_model import ElasticNetCV
from sklearn.datasets import make_regression

from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

############################
############################
# Function   #
############################
############################



##############
## data load

## for CH data
def load_changhai_drug_sen(ch_drug_sen_path):
    # ch_drug_sen_path = '/y/Bondi/data/changhai_data/PDC体外用药数据汇总：持续更新ing20200722(1).xlsx'
    ch_drug_sen_df = pd.read_excel(ch_drug_sen_path, header=2, index_col=0).iloc[:,3:]
    return ch_drug_sen_df

## for GDSC CCLE data
def get_gdsc_drug_sen(gdsc_bulk_data,drug=None,tissue=None,sp_sample_ls=None,drug_sen_type='ic50'):
    '''
    return Series: CELL_LINE_NAME, drug_sen
    '''
    if type(gdsc_bulk_data) == type('string'):
        gdsc_drug_sen_df = utilities.file_tools.read_txt(
            gdsc_drug_sen_path, columns_name_line=0)
    elif type(gdsc_bulk_data) == pd.core.frame.DataFrame:
        gdsc_drug_sen_df = gdsc_bulk_data
    else:
        return 'Wrong input gdsc_bulk_data'
    
    if drug is not None:
        sp_drug_index = gdsc_drug_sen_df.loc[gdsc_drug_sen_df['DRUG_NAME']==drug].index
    if tissue is not None:
        sp_tissue_index = gdsc_drug_sen_df.loc[gdsc_drug_sen_df['TCGA_DESC']==tissue].index
        
    select_index = set(sp_drug_index) & set(sp_tissue_index)
    
    if sp_sample_ls is not None:
        select_index = set(select_index) & set(sp_sample_ls)
    
    out_Series = None
    if drug_sen_type == 'ic50':
        out_Series = gdsc_drug_sen_df.loc[select_index].set_index('CELL_LINE_NAME')['LN_IC50']
    elif drug_sen_type == 'auc':
        out_Series = gdsc_drug_sen_df.loc[select_index].set_index('CELL_LINE_NAME')['AUC']
    
    return out_Series.sort_values().astype('float')

def load_GDSC_CCLE_data():
    pass

## for GTEx
def load_GTEx_gct_file(gct_path, dup_kep='first'):
    %R -i gct_path
    %R rt<-read.table(gct_path, skip = 2, header = TRUE, sep = "\t",check.names=F)
    GTEx_exp_mat = %R rt
    GTEx_exp_mat = GTEx_exp_mat.set_index('Description')
    GTEx_exp_mat = GTEx_exp_mat.drop('Name',axis=1)
    GTEx_exp_mat = GTEx_exp_mat[~GTEx_exp_mat.index.duplicated(keep=dup_kep)]
    return GTEx_exp_mat

#################
## data output
import zipfile
def make_zip(zip_filename, to_zip_ls):
    for fn in to_zip_ls:
        with zipfile.ZipFile(zip_filename, mode="a") as f:
            f.write(fn)          #追加写入压缩文件


#####################
## label transform
def XX():
    tpm_path = '/y/home/lgh/downloads/data/CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt'  ### <-
    tpm_exp_mat = pd.read_csv(tpm_path, sep='\t')
    %R library("AnnotationDbi")
    %R library("org.Hs.eg.db")
    ENSEMBL = np.array([i.split('.')[0] for i in tpm_exp_mat['gene_id']])
    %R -i ENSEMBL
    mapIds_ans = %R mapIds(org.Hs.eg.db,keys=ENSEMBL,column="SYMBOL",keytype="ENSEMBL",multiVals="first")

    new_tpm_exp_mat = tpm_exp_mat.loc[mapIds_ans !='NA']
    new_tpm_exp_mat.index = mapIds_ans[mapIds_ans !='NA']

    new_tpm_exp_mat = new_tpm_exp_mat.drop(['gene_id','transcript_id'],axis=1)
    new_tpm_exp_mat = new_tpm_exp_mat.loc[~new_tpm_exp_mat.index.duplicated(keep='first')]

##############
## data clean
def saple_name_transform(sample_name):
    return sample_name.replace('X','').replace('.','-')

def get_top_var_exp_mat():
    pass

def clean_drug_value(drug_value):
    try:
        return float(drug_value)
    except:
        return float(drug_value[1:])

def clean_data(expr_mat, dup='avg'):
    '''
    Remove duplicated row according to their index and columns content  
    Remove the abnormal sample or gene  
    input:  
    - annotation dataframe and expression matrix
    - dup = avg/plus/first/discard
    return: a cleaned expr_mat  
    '''
    # deal with duplicated
    no_duplicated_expr_mat = None
    duplicated_expr_mat = None

    # gene_list = df_annotation.loc[expr_mat.index,'gene']
    sample_list = expr_mat.columns

    dup_row = expr_mat.index.duplicated(keep=False)
    no_duplicated_expr_mat = expr_mat[~dup_row]
    duplicated_expr_mat = expr_mat[dup_row]
    
    # 对重复基因提供四种处理策略:平均值,相加,第一值,删掉
    if dup == 'avg':
        ## avg:
        duplicated_gene_set_ls = list(set(duplicated_expr_mat.index))
        avg_duplicated_expr_mat = np.zeros([len(duplicated_gene_set_ls), len(sample_list)])
        for i_g in range(len(duplicated_gene_set_ls)):
            for i_s in range(len(sample_list)):
                avg_duplicated_expr_mat[i_g][i_s] = np.mean(
                    duplicated_expr_mat.loc[duplicated_gene_set_ls[i_g]][sample_list[i_s]]
                    )
        
        avg_duplicated_expr_mat = pd.DataFrame(
            avg_duplicated_expr_mat, 
            index=duplicated_gene_set_ls, 
            columns=sample_list)
    
        unduplicated_expr_mat = pd.concat(
            [no_duplicated_expr_mat, avg_duplicated_expr_mat], axis=0)

    elif dup == 'plus':    
        ## plus
        pass
    
    elif dup == 'plus':  
        ## first
        pass

    elif dup == 'plus':  
        ## discard
        pass

    return expr_mat[~expr_mat.index.duplicated(keep='first')]

def filter_by_var(expr_mat, which='gene', what='nozero'):
    '''
    exp_mat: gene*sample
    which = 'gene'
    what = nozero or top_number
    '''
    if what != 'nozero': # top var
        expr_mat = expr_mat.loc[np.sum(expr_mat, axis=1) != 0]
        return expr_mat.loc[np.var(expr_mat,axis=1).sort_values().iloc[-what:].index]
    else: # var > 0        
        return expr_mat.loc[abs(np.var(expr_mat,axis=1))>0.000001]

def get_sample_with_sp_mut(maf_df,gene_symbol,ITDorNOT=None):
    if ITDorNOT is None:
        sp_maf_df = maf_df.loc[maf_df['Hugo_Symbol']==gene_symbol]
    else:
        sp_maf_df = maf_df.loc[maf_df['Hugo_Symbol']==gene_symbol].loc[Variants_df['ITDorNOT']=='ITD']
    res_ls = list(set(sp_maf_df['Tumor_Sample_Barcode']))
    print('number of sample:',len(res_ls))
    
    return res_ls

# only for TCGA
def get_type_of_TCGAsample(s_n):
    if int(s_n.split('-')[3][:2])<10:
        return 'Cancer'
    else:
        return 'Normal'

def TCGA_data_clean(exp_mat,id_symbol_anno):
    '''
    input：
    
    output：
    
    
    去掉'_’的行
    去掉全0行
    获取每个样本是病例还是对照
    '''
    exp_mat = exp_mat.loc[np.array(pd.Series(exp_mat.index).apply(lambda x:x[0]!='_'))]
    exp_mat = exp_mat.loc[np.sum(exp_mat,axis=1)!=0]
    gene_symbol = id_symbol_anno.loc[exp_mat.index]
    exp_mat.index = gene_symbol
    exp_mat = exp_mat[~exp_mat.index.duplicated(keep='first')]
    exp_mat = np.exp2(exp_mat)-1
    return {
        'exp_mat':exp_mat,
        'type_ls':list(map(get_type_of_TCGAsample,exp_mat.columns))
    }

############################
# Data Translation
def do_z_tran(exp_mat=None, series=None):
    '''
    exp_mat: gene*sample
    对每个基因进行的，使得每个基因在样本中的表达正太分布
    '''
    if exp_mat is not None:
        new_mat = preprocessing.scale(exp_mat,axis=1, with_mean=True, with_std=True, copy=True)
        new_mat = pd.DataFrame(new_mat,index=exp_mat.index, columns=exp_mat.columns)
        return new_mat
    elif series is not None:
        new_series = preprocessing.scale(series,axis=0, with_mean=True, with_std=True, copy=True)
        new_series = pd.Series(new_series,index=series.index)
        return new_series
    else:
        print('!!!')
        return '!!!'


############################
## stat
def do_cor(x, y, x_name=None, y_name=None, c=None, plot=False, figsize=None, plot_text=None):
    '''
    TO add draw line function
    '''
    if plot:
        fig, ax = plt.subplots(figsize=figsize)
        ax0 = ax.scatter(x, y, c=c)
        if c:
            plt.colorbar(ax0)

        ax.set_xlabel(x_name, fontsize=15)
        ax.set_ylabel(y_name, fontsize=15)
        ax.set_title('%s : %s~%s'%(plot_text, y_name, x_name))
        plt.show()

    %R -i x
    %R -i y
    %R lm_res = lm(x~y)
    %R lm_summary = summary(lm_res)

    coef = %R lm_summary$coefficients
    adj_r_squared = %R lm_summary$adj.r.squared
    r_squared = %R lm_summary$r.squared
    pvalue = coef[-1][-1]
    adj_r_squared = adj_r_squared[0]
    r_squared = r_squared[0]

    if plot:
        %R print(lm_summary)

    return {
        'pvalue':pvalue,
        'adj_r_squared':adj_r_squared,
        'r_squared':r_squared,
    }

def do_multi_lr(X_df,y):
    '''
    X_df: DataFrame 每列是每个feature的数值
    y: np.array 或者 Series
    '''
    # train_x = do_z_tran(exp_mat_deseq2_norm.loc[fry_feature_ls][sen_ls+unsen_ls]).T
    # train_y = do_z_tran(series=drug_sen_series[sen_ls+unsen_ls])
    %R -i X_df
    %R -i y
    %R lm_res <- lm(y~.,data = X_df)
    %R print(summary(lm_res))
    coef = %R lm_res$coefficients # intercept + other feature
    lm_res = %R lm_res
    return {'coef':coef,'lm_res':lm_res}

############################
## plot
def plot_bar_scatter(df_for_plot):
    '''
    df_for_plot = pd.DataFrame({
        'gene_ls_for_plot':gene_ls_for_plot,
        'avg_drug_sen_ls':avg_drug_sen_ls,
        'std_err_ls':std_err_ls,
        'ls_for_plot':ls_for_plot}).set_index('gene_ls_for_plot').sort_values('avg_drug_sen_ls')
    
    ls_for_plot is the list of list of sample for every gene to plot
    '''
    n = df_for_plot.shape[0]

    ### bar plot
    gene_ls_for_plot = list(df_for_plot.index)
    x = np.arange(n)
    y = list(df_for_plot['avg_drug_sen_ls'])
    std_err = list(df_for_plot['std_err_ls'])

    font2 = {
            'weight' : 'normal',
            'size' : 15,
            }

    error_params=dict(elinewidth=2,ecolor='black',capsize=3)#设置误差标记参数
    #绘制柱状图，设置误差标记以及柱状图标签

    fig, ax = plt.subplots(figsize=[12,5])

    ax.bar(x,y,yerr=std_err,error_kw=error_params, zorder=0)
    ax.set_xticks(np.arange(len(gene_ls_for_plot)))
    ax.set_xticklabels(gene_ls_for_plot)
    ax.set_xlabel('common mutation',font2)
    ax.set_ylabel('AUC',font2)

    ### scatter
    scatter_x = []
    scatter_y = []
    for i_d in range(0,n): ###################
        for c in df_for_plot['ls_for_plot'].iloc[i_d]:
            #x = np.random.normal(i, 0.02)
            scatter_x.append(np.random.normal(i_d, 0.1))
            #for_plot_x.append(i_d)
            scatter_y.append(c)

    ax.scatter(scatter_x, scatter_y,marker='*',color='r',s=1, zorder=1)

    ### post plot
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
            rotation_mode="anchor")
    #显示图形
    ax.set_title('%s'%(text))
    plt.show()


# wrap pheatmap
def plot_pheatmap(exp_mat, value_series=None, scale='none'):
    '''
    exp_mat : gene*sample

    exp_mat = exp_mat_deseq2_norm[sen_ls+unsen_ls].loc[tmp_feature_ls]
    scale: columns/none

    * 注意pheatmap的参数设置

    return: 
        cluster_array 是一个array，与exp.columns同长，每一个元素的值代表此saample被分到哪一个cluster
    '''
    mat_for_plot = exp_mat.T
    %R library(pheatmap)
    %R -i scale
    %R -i mat_for_plot
    %R p_res = pheatmap(mat_for_plot,scale=scale,cutree_row = 2,clustering_method = "ward.D",clustering_distance_rows = "correlation")
    order = %R p_res$tree_row$order
    labels = %R p_res$tree_row$label
    
    if value_series is not None:
        plt.scatter(
            np.arange(len(labels)),
            value_series.loc[mat_for_plot.index].loc[labels[order-1]],
                )

    cluster_array = %R cutree(p_res$tree_row,k=2)

    return cluster_array


## wrap ggolpt2 ##
def ggplot2_boxplot(df, x_level, y_value, stat_method=np.nan, ref_group = ".all.", title=None, pairwise=False):
    '''
    input:
        df: 2 columns: value, level
            ! 注意，level列应该为factor

        stat_method: "t.test" or np.nan or other...
    '''
    if title is None:
        title = 'Plot: %s ~ %s'%(y_value,x_level)
    
    %R library(ggplot2)
    %R library(ggpubr)        
    %R -i df
    %R -i x_level
    %R -i y_value
    %R -i stat_method
    %R if (is.nan(stat_method)){stat_method=NULL}
    %R -i ref_group
    %R -i title
    %R p <- ggboxplot(df, x = x_level, y = y_value,color = x_level, add = "jitter")
    %R p = p+labs(title=title)
    #if stat_method is not np.nan:
    if pairwise:
        y_hline = np.mean(df[y_value])
        %R -i y_hline
        %R plot(p + stat_compare_means(method=stat_method,label = "p.signif",ref.group=ref_group,hide.ns = FALSE)+geom_hline(yintercept = y_hline, linetype = 2))
        # %R plot(p + stat_compare_means(method=stat_method,label = "p.signif",ref.group=ref_group,hide.ns = FALSE)
        #%R plot(p+geom_hline(yintercept = y_hline, linetype = 2))
    else:
        %R plot(p + stat_compare_means(stat_method))
    #else:
    #    %R plot(p + stat_compare_means())
    print('ggplot2!')

def ggplot_violin(df,x_level, y_value):
    %R e <- ggplot(ToothGrowth, aes(x = dose, y = len))
    %R e = e + geom_violin()
    %R e = e + geom_violin(trim = FALSE,aes(fill = dose)) + geom_boxplot(width = 0.2)
    %R plot(e)
    print('')    

def ggplot2_volcano(deseq2_res, th_p=0.05, th_lfc=1):
    '''
    deseq2_res: columns: [log2FoldChange, padj]
    '''
    %R library("labeling")
    %R library(ggplot2)
    %R -i deseq2_res
    %R -i th_p
    %R -i th_lfc
    #将符合要求的筛出来
    %R threshold<-as.factor((deseq2_res$log2FoldChange>th_lfc|deseq2_res$log2FoldChange<(-th_lfc)) & deseq2_res$padj<th_p )
    %R p = ggplot(deseq2_res,aes(x= log2FoldChange,y= (-1)*log10(deseq2_res$padj),colour=threshold))+xlab("log2 fold-change")+ylab("-log10 p-value")+geom_point()
    %R plot(p)
    print('ggplot2_volcano')



# ↑ ↑ ↑ 

####################
# Genomic analysis
def cal_tmb(maf_path,captureSize=50,logScale = TRUE):
    '''
    %R -i maf_df
    %R -i captureSize
    %R -i logScale
    tmb_df = %R tmb(maf = maf_df,captureSize = captureSize,logScale = logScale)
    tmb_df = tmb_df.set_index('Tumor_Sample_Barcode')
    return tmb_df
    '''
    %R library(tidyverse)
    %R library(maftools)
    %R -i maf_path
    %R -i captureSize
    if logScale:
        %R logScale=TRUE
    else:
        %R logScale=FALSE

    %R library(maftools)
    %R library(tidyverse)
    %R rt <- read.maf(maf_path)
    %R maf = tmb(maf = rt, captureSize = 50, logScale = TRUE)
    return tmb_df

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


####################
# RNA analysis

### WGCNA
def do_sample_name_transform(moduleEigengenes_df):
    sample_ls = list(moduleEigengenes_df.index)
    sample_ls = list(map(saple_name_transform, sample_ls))
    moduleEigengenes_df.index = sample_ls
    return moduleEigengenes_df

def do_WGCNA(expr_mat, n_top_mad=2000):
    '''
    # use rpy2.robjects instead of %R
    import:
        expr_mat
    output:
        sft,moduleColors,colors,moduleEigengenes_df
    '''
    """    
    %R library(WGCNA)
    %R -i expr_mat
    %R WGCNA_matrix = t(expr_mat[order(apply(expr_mat, 1, mad), decreasing = T)[1:2000],])
    WGCNA_gene_ls = %R colnames(WGCNA_matrix)

    %R powers = c(c(1:10), seq(from = 12, to=20, by=2))
    %R sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 1)

    %R net = blockwiseModules(WGCNA_matrix,power = sft$powerEstimate,maxBlockSize = 6000,TOMType = "unsigned",minModuleSize = 30,reassignThreshold = 0,mergeCutHeight = 0.25,numericLabels = TRUE,pamRespectsDendro = FALSE,saveTOMs = F,verbose = 1)
    %R print('blockwiseModules DONE!!')
    %R moduleColors <- labels2colors(net$colors)
    %R MEs0 = moduleEigengenes(WGCNA_matrix, moduleColors)$eigengenes
    %R MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩 (样本vs模块)

    moduleColors = %R moduleColors
    colors = %R net$colors
    moduleEigengenes_df = %R MEs
    """
    pandas2ri.activate() ## 对于rpy2转pandas很有用

    with localconverter(ro.default_converter + pandas2ri.converter):
        r_expr_mat = ro.conversion.py2rpy(expr_mat)

    rcode = """
    expr_mat <- %s
    library(WGCNA)
    WGCNA_matrix = t(expr_mat[order(apply(expr_mat, 1, mad), decreasing = T)[1:2000],])

    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 1)

    net = blockwiseModules(WGCNA_matrix,power = sft$powerEstimate,maxBlockSize = 6000,TOMType = "unsigned",minModuleSize = 30,reassignThreshold = 0,mergeCutHeight = 0.25,numericLabels = TRUE,pamRespectsDendro = FALSE,saveTOMs = F,verbose = 1)
    print('blockwiseModules DONE!!')
    moduleColors <- labels2colors(net$colors)
    MEs0 = moduleEigengenes(WGCNA_matrix, moduleColors)$eigengenes
    MEs = orderMEs(MEs0)
    """%(r_expr_mat.r_repr())

    ro.r(rcode)

    WGCNA_gene_ls = ro.r("colnames(WGCNA_matrix)")
    sft = ro.r("sft")
    moduleColors = ro.r("moduleColors")
    colors = ro.r("net$colors")
    moduleEigengenes_df = ro.r("MEs")


    return {
            'WGCNA_gene_ls':WGCNA_gene_ls,
            'sft':sft,
            'moduleColors':moduleColors,
            'colors':colors,
            'moduleEigengenes_df':moduleEigengenes_df
            }

def WGCNA_analysis(expr_mat, drug_sen_series, n_top_mad=2000, plot_flag=True, samp_nam_trans=False):
    '''
    data:
    expr_mat:
    drug_sen_series:

    func:
    do_sample_name_transform()
    do_cor()
    '''

    do_WGCNA_res = do_WGCNA(expr_mat, n_top_mad)

    colors = do_WGCNA_res['colors']
    moduleColors = do_WGCNA_res['moduleColors']
    moduleEigengenes_df = do_WGCNA_res['moduleEigengenes_df']
    WGCNA_gene_ls = do_WGCNA_res['WGCNA_gene_ls']

    cluster_df = {'gene':[],'moduleColors':[],'color_id':[]}
    for i in range(len(colors)):
        cluster_df['gene'].append(WGCNA_gene_ls[i])
        cluster_df['moduleColors'].append(moduleColors[i])
        cluster_df['color_id'].append(colors[i])
    cluster_df = pd.DataFrame(cluster_df).set_index(['moduleColors'])

    n_cluster = len(set(colors))

    clusters = [[] for i in range(n_cluster)]
    for i in range(len(colors)):
        clusters[int(colors[i])].append(WGCNA_gene_ls[i])
    clusters = dict([[i, clusters[i]] for i in range(len(clusters))])

    ## get moduleEigengenes of each cluster
    if samp_nam_trans:
        moduleEigengenes_df = do_sample_name_transform(moduleEigengenes_df)

    ## do corr
    auc = np.array(drug_sen_series.loc[moduleEigengenes_df.index])

    moduleEigengenes_drug_cor_dict = {'moduleEigengenes':[],'Adjusted R-squared':[],'Multiple R-squared':[],'p-value':[]}
    for m in moduleEigengenes_df.columns:
        expr = np.array(moduleEigengenes_df[m])
        
        do_cor_res = do_cor(expr, auc, x_name='%s (rpkm)'%m, y_name='AUC', plot=True, plot_text='scatter')

        '''    
        plt.scatter(expr, auc)
        plt.title('%s'%m)
        
        %R -i expr
        %R -i auc

        corvalue = %R cor(expr, auc)
        %R lm_res = lm(auc~expr)
        plt.show()
        %R print(summary(lm_res))  
        
        %R summary_lm_res = summary(lm_res)
        adj_r_squared = %R summary_lm_res$adj.r.squared
        r_squared = %R summary_lm_res$r.squared
        coef = %R summary(lm_res)$coefficients
        p_value = coef[1][-1]
        '''
        moduleEigengenes_drug_cor_dict['moduleEigengenes'].append(m)
        moduleEigengenes_drug_cor_dict['Adjusted R-squared'].append(do_cor_res['adj_r_squared'])
        moduleEigengenes_drug_cor_dict['Multiple R-squared'].append(do_cor_res['r_squared'])
        moduleEigengenes_drug_cor_dict['p-value'].append(do_cor_res['pvalue'])

    moduleEigengenes_drug_cor_df = pd.DataFrame(moduleEigengenes_drug_cor_dict).set_index('moduleEigengenes')

    return {
            'moduleEigengenes_drug_cor_df':moduleEigengenes_drug_cor_df,
            'do_WGCNA_res':do_WGCNA_res,
            'cluster_df':cluster_df
        }
## WGCNA ## ↑

### enrichment analysis
def do_ssGSEA():
    pass

def do_GSVA():
    pass

def do_GSEA():
    pass

def do_enrichplot(feature_array):
    '''
    GOcir_enrich_result_df
    enrich_result_df : for barplot() : use R script
        ```
        go <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.2,keyType = 'ENTREZID')
        barplot(go,showCategory=20,drop=T)
        ```
    '''
    %R library("clusterProfiler")
    %R library("org.Hs.eg.db")
    %R library("enrichplot")
    %R library("ggplot2")

    %R -i feature_array
    # %R DEG_gene_array <- ry_feature_array
    %R df2 <- bitr(feature_array,fromType = "SYMBOL",toType = c("ENTREZID", "GENENAME"),OrgDb = org.Hs.eg.db)
    %R feature_gene_eid = df2[[2]]
    %R kk=enrichGO(gene=feature_gene_eid,OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, ont="ALL", readable=T)
    enrich_result_df = %R as.data.frame(kk)
    '''
    GOcir_enrich_result_df = enrich_result_df[['ONTOLOGY','ID','Description','geneID','p.adjust']]
    GOcir_enrich_result_df.columns = ['Category','ID','Term','Genes','adj_pval']
    GOcir_enrich_result_df['Genes'] = enrich_result_df['Genes'].map(lambda x: x.replace('/',','))
    GOcir_enrich_result_df.index = list(range(1,len(enrich_result_df)+1))   ###
    '''


    return {'GOcir_enrich_result_df':GOcir_enrich_result_df,'enrich_result_df':enrich_result_df}

def prepare_for_GOcircle(feature_sorted_df ,enrich_result_df, method, save_path):
    '''
    feature_sorted_df: two columns: first is gene, second is the value to be sort (maybe the logFC)
    enrich_result_df: result of do_enrichplot
    '''
    feature_sorted_df.columns = ['ID','logFC']

    feature_sorted_df.to_csv('%s/%s_for_GOPlot.csv'%(save_path,method),sep='\t')
    enrich_result_df.to_csv('%s/enrich_result_for_GOcircle_%s_df.tsv'%(save_path,method),sep='\t')


############################
## feature selection
def do_deseq2(expr_mat, condition_table, colData):
    '''
    expr_mat : genes*samples
    '''
    %R -i expr_mat
    %R -i condition_table
    %R as.factor(condition_table)

    %R library(DESeq2)
    %R dds <- DESeqDataSetFromMatrix(expr_mat, DataFrame(condition_table), design= ~ condition_table )
    %R dds2 <- DESeq(dds)
    %R res <- results(dds2)
    DESeq_res = %R data.frame(res)
    %R rld <- vst(dds2, blind = T)
    %R expr_norm <- assay(rld)
    expr_norm = %R expr_norm

    DESeq_res.index = expr_mat.index
    expr_norm_mat = pd.DataFrame(expr_norm, index=expr_mat.index, columns=expr_mat.columns)

    return {'DESeq_res':DESeq_res, 'expr_norm_mat':expr_norm_mat}

def cor_onebyone(expr_mat, drug_sen_array):
    '''
    expr_mat : genes × samples ， maybe rpkm
    drug_sen_array ：length = samples
    '''
    corr_res_dict = {'gene':[], 'adj_r_squared':[], 'pvalue':[]}

    for gene in expr_mat.index:
        if type(expr_mat.loc[gene]) is pd.core.series.Series:
            expr = np.array(expr_mat.loc[gene])
        else:
            expr = np.array(expr_mat.loc[gene].iloc[0])
        
        %R -i expr
        %R -i drug_sen_array
        #cor_value = %R cor(expr, auc)
        #cor_value = cor_value[0]
        %R lm_res = lm(drug_sen_array ~ expr)
        %R lm_summary = summary(lm_res)
        coef = %R lm_summary$coefficients
        adj_r_squared = %R lm_summary$adj.r.squared
        pvalue = coef[-1][-1]
        adj_r_squared = adj_r_squared[0]
        
        corr_res_dict['gene'].append(gene)
        corr_res_dict['adj_r_squared'].append(adj_r_squared)
        corr_res_dict['pvalue'].append(pvalue)
        #if adj_r_squared>0.25 and pvalue<0.05:
        #    gene_ls.append(gene)

    corr_res_df = pd.DataFrame(corr_res_dict).set_index(['gene'])

    return corr_res_df

def do_elastic_net(X, y, gene_ls, alphas, cv, max_iter, random_state, selection, verbose=1):
    '''
    X: sample*gene np.array
    y: sample np.array
    gene_ls:np.array
    '''
    regr = ElasticNetCV(alphas=alphas, cv=cv, max_iter=max_iter, random_state=random_state, selection=selection)
    # print(X.shape)
    regr.fit(X, y)

    # predict_y = regr.predict(X)
    '''
    feature_ls = []
    feature_weight_df = []
    for i in range(len(regr.coef_)):
        if regr.coef_[i]!=0:
            feature_ls.append(all_feature_ls[i])
            feature_weight_df.append([all_feature_ls[i], regr.coef_[i]])
    feature_weight_df = pd.DataFrame(feature_weight_df,columns=['Feature','Weight']).set_index(['Feature'])
    '''
    feature_weight_df = pd.DataFrame({'Feature':gene_ls[regr.coef_!=0],
                                    'Weight':regr.coef_[regr.coef_!=0]}).set_index(['Feature'])
    feature_ls = np.array(feature_weight_df.index)


    if verbose:  
        print('Elastic Net keep %s features...'%len(feature_ls))
    
    return {'model':regr, 'feature_ls':feature_ls, 
            'feature_weight_df':feature_weight_df, 
            'X':X, 'y':y}


from sklearn.ensemble import RandomForestRegressor
def do_random_forest(X, y, all_feature_ls, max_depth, random_state, verbose=1):
    '''
    X
    y
    all_feature_ls

    max_depth
    random_state
    '''
    regr = RandomForestRegressor(max_depth=max_depth, random_state=random_state)
    regr.fit(X, y)

    predict_y = regr.predict(X)

    feature_ls = []
    feature_weight_df = []
    for i in range(len(regr.feature_importances_)):
        if regr.feature_importances_[i]!=0:
            feature_ls.append(all_feature_ls[i])
            feature_weight_df.append([all_feature_ls[i], regr.feature_importances_[i]])
    feature_weight_df = pd.DataFrame(feature_weight_df,columns=['Feature','Weight']).set_index(['Feature'])
    
    if verbose:  
        print('Random Forest keep %s features...'%len(feature_ls))
    
    return {'model':regr, 'feature_ls':feature_ls, 
            'feature_weight_df':feature_weight_df, 
            'X':X, 'y':y}

### pheatmap 和 逐步选择
def step_wise_select_by_pheatmap(exp_mat, condition_table, candi_len=500):
    '''
    exp_mat: gene*sample
    condition_table: as long as exp_mat.columns
    
    exp_mat = exp_mat_deseq2_norm[sen_ls+unsen_ls].loc[DEG_gene_ls]
    condition_table = [0]*len(sen_ls)+[1]*len(unsen_ls)
    '''
    # DEG_gene_ls = do_deseq2_df.sort_values('padj').iloc[:500].index
    # condition_table = [0]*len(sen_ls)+[1]*len(unsen_ls)
    candi_feature_ls = list(exp_mat.index[:candi_len])

    tmp_feature_ls = [candi_feature_ls[0]]
    mi = 0
    for g in candi_feature_ls[1:len(candi_feature_ls)]:
        tmp_feature_ls.append(g)

        # mat_for_plot = res_do_deseq2['expr_norm_mat'].loc[tmp_feature_ls].T
        mat_for_plot = exp_mat.loc[tmp_feature_ls].T
        %R library(pheatmap)
        %R -i mat_for_plot
        %R p_res = pheatmap(mat_for_plot,scale='column',cutree_row = 2,clustering_method = "ward.D",clustering_distance_rows = "correlation",silent=TRUE)
        # order = %R p_res$tree_row$order
        # labels = %R p_res$tree_row$label
        row_cluster = %R cutree(p_res$tree_row,k=2)
        tmp_mi = metrics.adjusted_mutual_info_score(condition_table, row_cluster)
        #print('mi now is %s'%tmp_mi)
        if tmp_mi < mi:
            #print('drop')
            tmp_feature_ls = tmp_feature_ls[:-1]
        else:
            mi = tmp_mi

    print(mi)

    selected_feature_ls = tmp_feature_ls
    return selected_feature_ls

def plot_pheatmap_drug_sen(exp_mat, drug_sen_series, condition_table, scale='none'):
    '''
    exp_mat : gene*sample
    condition_table : the same as len of sample

    exp_mat = exp_mat_deseq2_norm[sen_ls+unsen_ls].loc[tmp_feature_ls]
    scale: columns/none
    '''
    mat_for_plot = exp_mat.T
    %R library(pheatmap)
    %R -i scale
    %R -i mat_for_plot
    %R p_res = pheatmap(mat_for_plot,scale=scale,cutree_row = 2,clustering_method = "ward.D",clustering_distance_rows = "correlation")
    order = %R p_res$tree_row$order
    labels = %R p_res$tree_row$label
    plt.scatter(
        np.arange(len(labels)),
        drug_sen_series.loc[mat_for_plot.index].loc[labels[order-1]],
                )

    row_cluster = %R cutree(p_res$tree_row,k=2)
    mi = metrics.adjusted_mutual_info_score(condition_table, row_cluster)

    #return {'group_df':group_df,
    #        'mi':mi
    #        }




## 整合的 pipleline function

def do_diff_corr_wgcna(drug, all_drug_sen_df, rpkm_expr_mat, get_sen_unsen_p , get_sen_unsen_c):
    '''
    drug
    all_drug_sen_df
    rpkm_expr_mat
    get_sen_unsen_p
    get_sen_unsen_c

    # example: /y/Bondi/jupyter/hl/4-1-3分析/for_hl.ipynb
    '''

    '''
    i_drug = 3  ### can be changed

    print(drug[i_drug])

    # get data
    
    drug_sen_path = '/y/Bondi/data/VIZOME/nature_aml_drug_sen.tsv'
    all_drug_sen_df = pd.read_csv(drug_sen_path, sep='\t', index_col=0)

    rpkm_path = '/y/Bondi/data/nature_aml_fpkm.txt'
    rpkm_expr_mat = pd.read_csv(rpkm_path, sep='\t', index_col=0)
    '''
    ## get specific drug sen
    drug_sen_df = all_drug_sen_df.loc[drug].set_index(['lab_id']) # can change
    # get useful
    useful_cell = set(rpkm_expr_mat.columns) & set(drug_sen_df.index) # can change

    # filter sample by useful, filter gene by is duplicate
    rpkm_expr_mat = rpkm_expr_mat[useful_cell]
    rpkm_expr_mat = clean_data(rpkm_expr_mat, dup='avg')
    drug_sen_df = drug_sen_df.loc[useful_cell]

    # 获取上下20% : get sen unsen
    get_sen_unsen_res = get_sensitivity_unsensitivity(drug_sen_df, useful_cell, p=0.2, c=None, indicator='auc')
    diff_cell = get_sen_unsen_res['sensitivity_cell_ls']+get_sen_unsen_res['unsensitivity_cell_ls']
    diff_rpkm_expr_mat = rpkm_expr_mat[get_sen_unsen_res['sensitivity_cell_ls']+get_sen_unsen_res['unsensitivity_cell_ls']]

    # wgcna

    WGCNA_analysis_res = WGCNA_analysis(rpkm_expr_mat, drug_sen_df, n_top_mad=2000, plot_flag=True)

    # limma
    condition_table = np.array(['sen']*len(get_sen_unsen_res['sensitivity_cell_ls'])+
                        ['unsen']*len(get_sen_unsen_res['unsensitivity_cell_ls']))
    do_limma_res = do_limma(expr_mat=diff_rpkm_expr_mat, 
            condition_table=condition_table, 
            logFC_th=1, pvalue_th=0.05)


    # cor
    cor_res_df = cor_onebyone(expr_mat=rpkm_expr_mat, 
                drug_sen_array=np.array(drug_sen_df['auc']))



    # interaction

    # len(set(WGCNA_signi_ls) & set(significent_gene_ls) & set(get_diff_from_limma_res['diff_gene_ls']))

    return {
        'WGCNA_analysis_res':WGCNA_analysis_res,
        'do_limma_res':do_limma_res,
        'cor_res_df':cor_res_df
    }


#####################################
# prepare for other software

## for GSEA

def prepare_exp_mat(expr_mat, condition_table, path, filename=):
    '''
    expr_mat: gene × sample
    path: path to save, tab分隔, .gct!
    '''
    expr_mat_forGSEA = pd.concat([expr_mat['Symbol'],expr_mat],axis=1)
    colnames = list(exp_mat_forGSEA.columns)
    colnames[0] = 'NAME'
    colnames[1] = 'Description'
    exp_mat_forGSEA.columns = colnames
    exp_mat_forGSEA.to_csv('%s/%s.gct'%(path,filename), sep='\t')
    # ===

    pheno_1 = condition_table[0]
    pheno_2 = (set(condition_table)-set([pheno_1])).pop()
    sample_n = len(condition_table)
    group_n = 2

    cls_file = ''
    cls_file += '%s\t%s\t1\n'%(sample_n, group_n)
    cls_file += '#\t%s\t%s\n'%(pheno_1,pheno_2)
    cls_file += '\t'.join(condition_table)

    f = open('%s/%s.cls'%(path,filename), 'w', encoding='utf-8') 
    f.write(cls_file)
    f.close()

    return exp_mat_forGSEA



#####################################
# 