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

import sys
sys.path.append("/y/Bondi/src_sync/")

import utilities.file_tools

def load_changhai_drug_sen(ch_drug_sen_path):
    # ch_drug_sen_path = '/y/Bondi/data/changhai_data/PDC体外用药数据汇总：持续更新ing20200722(1).xlsx'
    ch_drug_sen_df = pd.read_excel(ch_drug_sen_path, header=2, index_col=0).iloc[:,3:]
    return ch_drug_sen_df

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


def saple_name_transform(sample_name):
    return sample_name.replace('X','').replace('.','-')

def do_sample_name_transform(moduleEigengenes_df):
    sample_ls = list(moduleEigengenes_df.index)
    sample_ls = list(map(saple_name_transform, sample_ls))
    moduleEigengenes_df.index = sample_ls
    return moduleEigengenes_df

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

def clean_drug_df():
    pass

def clean_drug_value(drug_value):
    try:
        return float(drug_value)
    except:
        return float(drug_value[1:])

###################
## wrap function ##
###################
def data_filter(which, excluded_ls=None, included_ls=None):
    '''
    Input
        which: filter object is gene or sample
        excluded_ls: list or np.array 正则表达式??
        included_ls: list of np.array

    '''

def import_data(genomic_alter_path, expr_mat_path, drug_response_path):
    drug_response_df = pd.read_csv(drug_response_path,index_col=0)
    expr_mat = pd.read_csv(expr_mat_path,index_col=0)

    # remove deficiency data
    sample_list = list(set(drug_response_df.index) & set(expr_mat.columns))
    gene_list = list(expr_mat.columns)

    expr_mat = expr_mat.T.loc[sample_list]
    drug_response_df = drug_response_df.loc[sample_list]

    return {'expr_mat':expr_mat,
            'drug_response_df':drug_response_df,
            'sample_list':sample_list,
            'gene_list':gene_list
        }

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



def do_limma(expr_mat, condition_table, logFC_th=1, pvalue_th=0.05):
    '''
    expr_mat: genes × samples : rpkm
    condition_table: samples × 2
    '''
    %R library(edgeR)
    %R library(limma)

    #condition_table = np.array(['sen']*len(get_sen_unsen_res['sensitivity_cell_ls'])+
    #                    ['unsen']*len(get_sen_unsen_res['unsensitivity_cell_ls']))
    #log2_diff_expr_mat = np.log2(diff_expr_mat+0.1)

    %R -i expr_mat
    %R -i condition_table
    %R design <- model.matrix(~condition_table)
    %R colnames(design) <- levels(condition_table)
    %R rownames(design) <- colnames(expr_mat)
    %R fit <- lmFit(expr_mat, design)
    %R fit <- eBayes(fit, trend=TRUE)
    %R output <- topTable(fit, coef=2,n=Inf)
    # %R sum(output$adj.P.Val<0.05)
    limma_output = %R output

    gene_ls = np.array(limma_output.index)
    up_genes = list(set(gene_ls[limma_output['logFC']>logFC_th]) & set(gene_ls[limma_output['P.Value'] < pvalue_th]))
    down_genes = list(set(gene_ls[limma_output['logFC']<-logFC_th]) & set(gene_ls[limma_output['P.Value'] < pvalue_th]))

    return {
        'up_genes':up_genes,
        'down_genes':down_genes,
        'limma_output':limma_output
    }

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

def plot_pheatmap_drug_sen(exp_mat, condition_table):
    '''
    exp_mat : gene*sample
    condition_table : the same as len of sample

    exp_mat = exp_mat_deseq2_norm[sen_ls+unsen_ls].loc[tmp_feature_ls]
    '''
    mat_for_plot = exp_mat.T
    %R library(pheatmap)
    %R -i mat_for_plot
    %R p_res = pheatmap(mat_for_plot,scale='column',cutree_row = 2,clustering_method = "ward.D",clustering_distance_rows = "correlation")
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


# import exp data

'''
expr_mat_path = '/y/home/lgh/lgh/work/drug_sensitivity/data/all_breast_cell_com_geneexpression_norm.csv'
#drug_response_path = '/y/home/lgh/lgh/work/drug_sensitivity/data/ic50_formatted.csv'

import_data_res = import_data(genomic_alter_path=None,
        expr_mat_path=expr_mat_path,drug_response_path=drug_response_path)


expr_mat = import_data_res['expr_mat']
drug_response_df = import_data_res['drug_response_df']
'''

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



new_tpm_exp_mat.to_csv('/y/Bondi/data/CCLE/expression/CCLE_RNAseq_tpm.tsv',sep='\t')

# start from here

tpm_exp_mat = pd.read_csv('/y/Bondi/data/CCLE/expression/CCLE_RNAseq_tpm.tsv',sep='\t', index_col=0)
tpm_exp_mat.columns = [s.split('_')[0] for s in tpm_exp_mat.columns]
tpm_exp_mat = np.log(0.1+tpm_exp_mat)

## ccle annotation

## CCLE annotation
CCLE_annotation_file_path = '/y/Bondi/data/CCLE/sample_info/Cell_lines_annotations_20181226.txt'


CCLE_annotation_df = pd.read_csv(CCLE_annotation_file_path, sep='\t')
CCLE_annotation_df['CCLE_ID'] = [s.split('_')[0] for s in CCLE_annotation_df['CCLE_ID']]

NAME_CCLEID_annotation_df = CCLE_annotation_df.loc[CCLE_annotation_df['type_refined'] == CCLE_annotation_df['type_refined']].set_index('Name')
CCLEID_NAME_annotation_df = CCLE_annotation_df.loc[CCLE_annotation_df['type_refined'] == CCLE_annotation_df['type_refined']].set_index('CCLE_ID')


##########
# LFM-A13
##########
drug = 'LFM-A13'
drug_sen_path = '/y/Bondi/data/GDSC/cell_line_sen/LFM-A13.csv'

drug_sen_df = pd.read_csv(drug_sen_path).set_index('Cell line')
drug_sen_df = drug_sen_df.loc[set(NAME_CCLEID_annotation_df.index) & set(drug_sen_df.index)]
drug_sen_df.index = NAME_CCLEID_annotation_df.loc[drug_sen_df.index]['CCLE_ID']

drug_sen_series = drug_sen_df['IC50'].sort_values()
useful_sample = list(drug_sen_series.index)
sp_tpm_exp_mat = tpm_exp_mat[useful_sample]

# ==

## analysis
expr_mat = sp_tpm_exp_mat

sen_ls = list(drug_sen_series.iloc[:5].index)
unsen_ls = list(drug_sen_series.iloc[-5:].index)
condition_table=np.array([0]*5+[1]*5)

### 取top var

expr_mat = expr_mat.loc[np.sum(expr_mat, axis=1) != 0]
top_var_gene = list(np.var(expr_mat, axis=1).sort_values().iloc[-5000:].index)
expr_mat = expr_mat.loc[top_var_gene]

p_res = plot_pheatmap_drug_sen(exp_mat=expr_mat.loc[np.var(expr_mat, axis=1).sort_values().iloc[-100:].index][sen_ls+unsen_ls], 
                               condition_table=condition_table)

# ==

# DESeq2
do_limma_res = do_limma(expr_mat=expr_mat.loc[:][sen_ls+unsen_ls], 
         condition_table=condition_table, 
         logFC_th=1, pvalue_th=0.05)

DEG_ls = do_limma_res['limma_output'].sort_values('P.Value').iloc[:100].index
p_res = plot_pheatmap_drug_sen(exp_mat=expr_mat.loc[DEG_ls][sen_ls+unsen_ls], 
                               condition_table=condition_table)

do_limma_res['limma_output'].to_csv('./%s_limma_df.tsv'%drug,sep='\t')


# ==
#candi_gene_ls = do_limma_res['limma_output'].sort_values('P.Value').iloc[:500].index
candi_gene_ls = np.array( do_limma_res['limma_output'].loc[(abs(do_limma_res['limma_output']['logFC'])>1) & (do_limma_res['limma_output']['P.Value']<0.05)].index )

do_random_forest_res = do_random_forest(X=expr_mat.loc[candi_gene_ls].T, 
                 y=np.log(np.array(drug_sen_series)), 
                 all_feature_ls=np.array(expr_mat.loc[candi_gene_ls].index), 
                 max_depth=2, 
                 random_state=random.randint(100,200), 
                 verbose=1)


rf_gene_ls = np.array(list(do_random_forest_res['feature_weight_df'].index))

p_res = plot_pheatmap_drug_sen(exp_mat=expr_mat.loc[DEG_ls][sen_ls+unsen_ls], 
                               condition_table=condition_table)


# enrich
'''
deq2_DEG_df = pd.DataFrame({'ID':do_deseq2_df.loc[deq2_DEG]['log2FoldChange'].index,
              'logFC':list(do_deseq2_df.loc[deq2_DEG]['log2FoldChange'])})
deq2_DEG_df.to_csv('/y/Bondi/jupyter/hl/3-23-PY34/deq2_DEG_df_for_GOPlot.csv')
'''
%R -i rf_gene_ls
%R DEG_gene_array <- rf_gene_ls

%R df2 <- bitr(DEG_gene_array,fromType = "SYMBOL",toType = c("ENTREZID", "GENENAME"),OrgDb = org.Hs.eg.db)
DEG_gene_eid = %R df2[[2]]

%R -i DEG_gene_eid
%R kk=enrichGO(gene=DEG_gene_eid,OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, ont="ALL", readable=T)


enrich_result_df = %R as.data.frame(kk)
enrich_result_df = enrich_result_df[['ONTOLOGY','ID','Description','geneID','p.adjust']]
enrich_result_df.columns = ['Category','ID','Term','Genes','adj_pval']

enrich_result_df['Genes'] = enrich_result_df['Genes'].map(lambda x: x.replace('/',','))

enrich_result_df.index = list(range(1,len(enrich_result_df)+1))   ###

enrich_result_df.to_csv('./enrich_result_for_GOcircle_%s.tsv'%drug,sep='\t')
















