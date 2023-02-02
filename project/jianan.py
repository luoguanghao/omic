drug = ['GDC-0941', 'BEZ235', 'Venetoclax']


# get data
cpm_path = '/y/Bondi/data/VIZOME/nature_aml_count.tsv'
drug_sen_path = '/y/Bondi/data/VIZOME/nature_aml_drug_sen.tsv'
rpkm_path = '/y/Bondi/data/nature_aml_fpkm.txt'

drug_sen_df = pd.read_csv(drug_sen_path, sep='\t', index_col=0)
cpm_expr_mat = pd.read_csv(cpm_path, sep='\t', index_col=1).drop(['Gene'], axis=1)
rpkm_expr_mat = pd.read_csv(rpkm_path, sep='\t', index_col=1).drop(['Gene'], axis=1)


useful_cell = set(cpm_expr_mat.columns) & set(drug_sen_df.loc[drug[2]]['lab_id']) # can change

drug_sen_df = drug_sen_df.loc[drug[2]].set_index(['lab_id']) # can change

## 获取上下20%
get_sen_unsen_res = get_sensitivity_unsensitivity(drug_sen_df, useful_cell, p=0.2, c=None, indicator='ic50')
diff_cell = get_sen_unsen_res['sensitivity_cell_ls']+get_sen_unsen_res['unsensitivity_cell_ls']
diff_expr_mat = cpm_expr_mat[get_sen_unsen_res['sensitivity_cell_ls']+get_sen_unsen_res['unsensitivity_cell_ls']]

## get some result
auc_diff = np.array(drug_sen_2_df.loc[diff_cell]['auc'])

# diff exp
## try limma
%R library(edgeR)
%R library(limma)

condition_table = np.array(['sen']*len(get_sen_unsen_res['sensitivity_cell_ls'])+
                    ['unsen']*len(get_sen_unsen_res['unsensitivity_cell_ls']))
log2_diff_expr_mat = np.log2(diff_expr_mat+0.1)

%R -i log2_diff_expr_mat
%R -i condition_table
%R design <- model.matrix(~condition_table)
%R colnames(design) <- levels(condition_table)
%R rownames(design) <- colnames(log2_diff_expr_mat)
%R fit <- lmFit(log2_diff_expr_mat, design)
%R fit <- eBayes(fit, trend=TRUE)
%R output <- topTable(fit, coef=2,n=Inf)
%R sum(output$adj.P.Val<0.05)
limma_output = %R output

gene_ls = np.array(limma_output.index)
up_genes = list(set(gene_ls[limma_output['logFC']>np.log2(3)]) & set(gene_ls[limma_output['P.Value'] < 0.000001]))
down_genes = list(set(gene_ls[limma_output['logFC']<np.log2(0.3)]) & set(gene_ls[limma_output['P.Value'] < 0.000001]))

#################
# corr one by one

corr_res_dict = {'gene':[], 'adj_r_squared':[], 'pvalue':[]}

for gene in diff_expr_mat.index:
    if type(diff_expr_mat.loc[gene]) is pd.core.series.Series:
        expr = np.array(diff_expr_mat.loc[gene])
    else:
        expr = np.array(diff_expr_mat.loc[gene].iloc[0])
    
    %R -i expr
    %R -i auc
    #cor_value = %R cor(expr, auc)
    #cor_value = cor_value[0]
    %R lm_res = lm(auc~expr)
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

corr_signif_gene = []
for gene in corr_res_df.index:
    if corr_res_df.loc[gene]['adj_r_squared']>0.25 and corr_res_df.loc[gene]['pvalue']<5e-07:
        corr_signif_gene.append(gene)

for gene in set(corr_signif_gene) & set(up_genes+down_genes):
    print(gene)







#######################
# expression analysis #
#######################
import collections
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

%load_ext rpy2.ipython

import math, csv, pickle
import random, os


drug_ls = ['GDC-0941', 'BEZ235', 'PI-103', 'Venetoclax', 'Sorafenib']
i_drug = 3


def saple_name_transform(sample_name):
    return sample_name.replace('X','').replace('.','-')

def do_sample_name_transform(moduleEigengenes_df):
    sample_ls = list(moduleEigengenes_df.index)
    sample_ls = list(map(saple_name_transform, sample_ls))
    moduleEigengenes_df.index = sample_ls
    return moduleEigengenes_df

def old_get_sensitivity_unsensitivity(GDSC_df, cell_list, p=None, c=None, indicator='IC50'):
    '''
    取得药物敏感性排头和拍尾的细胞系名称
    对于药物敏感性，是IC50小的敏感性强， IC50大的敏感性不强
    
    percentage: 以前后多少比例作为敏感/不敏感 <0.5
    indicator: IC50 / AUC   
    '''
    sorted_drug_GDSC_df = GDSC_df.loc[cell_list].sort_values(by=indicator, axis=0) # 排序是从小到大的
    
    if p != None:
        n = len(cell_list)
        tmp_n = int(n * p)
    else:
        tmp_n = c
        
    sensitivity_cell_ls = list(sorted_drug_GDSC_df.index[:tmp_n])
    unsensitivity_cell_ls = list(sorted_drug_GDSC_df.index[-tmp_n:])
    
    return {
            'sensitivity_cell_ls': sensitivity_cell_ls,
           'unsensitivity_cell_ls':unsensitivity_cell_ls
           }


def get_sensitivity_unsensitivity(drug_sen_series, cell_list=None, p=None, c=None):
    '''
    取得药物敏感性排头和拍尾的细胞系名称
    对于药物敏感性，是IC50小的敏感性强， IC50大的敏感性不强
    
    percentage: 以前后多少比例作为敏感/不敏感 <0.5
    indicator: IC50 / AUC   
    '''
    if cell_list!= None:
        sorted_drug_sen_series = drug_sen_series.loc[cell_list].sort_values() # 排序是从小到大的
    else:
        sorted_drug_sen_series = drug_sen_series.sort_values() # 排序是从小到大的
        cell_list = list(sorted_drug_sen_series.index)
    
    if p != None:
        n = len(cell_list)
        tmp_n = int(n * p)
    else:
        tmp_n = c
        
    sensitivity_cell_ls = list(sorted_drug_sen_series.index[:tmp_n])
    unsensitivity_cell_ls = list(sorted_drug_sen_series.index[-tmp_n:])
    
    return {
            'senls': sensitivity_cell_ls,
           'unsen_ls':unsensitivity_cell_ls
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

def filter_var(mat=None, array=None, th=0.1):
    '''
    - 依据方差筛选，方差太小的不要
    - 输入mat，按照列做var，算出每列var，返回哪些True哪些False的一个array
        - 对于格式为gene×sample的 expr_mat，传进来前要先转置
    '''
    if array is not None:
        if np.var(array) < th:
            return False
        return True
    
    if mat is not None:
        ans = np.var(mat) < th
        return ans


def sort_by_list(sort_by_ls, to_sort_ls):
    sort_index = np.array(sort_by_ls).argsort()
    sorted_ls = [to_sort_ls[i] for i in sort_index]
    return sorted_ls


import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

from rpy2.robjects.conversion import localconverter

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


def get_diff_from_limma(limma_output_df, logFC_th=1, pvalue_th=0.05):
    up_gene_ls = []
    down_gene_ls = []
    for gene in limma_output_df.index:
        if limma_output_df.loc[gene, 'logFC']>logFC_th and limma_output_df.loc[gene, 'adj.P.Val']<pvalue_th:
            up_gene_ls.append(gene)
        if limma_output_df.loc[gene, 'logFC']<-logFC_th and limma_output_df.loc[gene, 'adj.P.Val']<pvalue_th:
            down_gene_ls.append(gene)
    diff_gene_ls = up_gene_ls + down_gene_ls

    print('  diff_gene_ls length is :',len(diff_gene_ls))

    return {
        'up_gene_ls':up_gene_ls,
        'down_gene_ls':down_gene_ls,
        'diff_gene_ls':diff_gene_ls
    }


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


def do_t_test(data_ls, condition_table, plot=False):
    %R -i condition_table
    %R -i data_ls
    %R res = t.test(data_ls ~ condition_table)
    pvalue = %R res$p.value
    statistic = %R res$statistic
    ttest_res = %R res

    if plot:
        pass

    return {'pvalue':pvalue,
            'statistic':statistic,
            'ttest_res':ttest_res}



def t_test_onebyone(expr_mat, sen_ls, unsen_ls):
    '''
    expr_mat : genes × samples ， maybe rpkm
    drug_sen_array ：length = samples
    condition_table : 定义每个数据的属性，分组
    '''
    ttest_res_dict = {'gene':[], 'pvalue':[]}

    for gene in expr_mat.index:
        data_ls = expr_mat.loc[gene, sen_ls+unsen_ls]
        condition_table = np.array(['sen']*len(sen_ls) + ['unsen']*len(unsen_ls))
        do_t_test_res = do_t_test(data_ls, condition_table)
        
        ttest_res_dict['gene'].append(gene)
        ttest_res_dict['pvalue'].append(do_t_test_res['pvalue'])
        '''        
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
        '''

    ttest_res_df = pd.DataFrame(ttest_res_dict).set_index(['gene'])

    return ttest_res_df   


def get_significient_from_cor(cor_output_df, r2_th=0.25, pvalue_th=0.05):
    pos_gene_ls = []
    neg_gene_ls = []
    significent_gene_ls = []
    for gene in cor_output_df.index:
        if cor_output_df.loc[gene, 'adj_r_squared']>=r2_th and cor_output_df.loc[gene, 'pvalue']<=pvalue_th:
            significent_gene_ls.append(gene)
    
    print('   %s significent_gen'%len(significent_gene_ls))
    return {
        'significent_gene_ls':significent_gene_ls
    }

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


def do_diff_corr_wgcna(drug, all_drug_sen_df, rpkm_expr_mat, get_sen_unsen_p , get_sen_unsen_c):
    '''
    import data:
        drug
        all_drug_sen_df
        rpkm_expr_mat
        get_sen_unsen_p
        get_sen_unsen_c
    function:
        ...
    output:
        ...
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
    ## ...
    print("""
    =============
    || prepare ||
    =============
    """)

    drug_sen_df = all_drug_sen_df.loc[drug].set_index(['lab_id']) # can change

    useful_cell = set(rpkm_expr_mat.columns) & set(drug_sen_df.index) # can change


    rpkm_expr_mat = rpkm_expr_mat[useful_cell]
    rpkm_expr_mat = clean_data(rpkm_expr_mat, dup='avg')

    drug_sen_df = drug_sen_df.loc[useful_cell]

    ## 获取上下20%
    get_sen_unsen_res = get_sensitivity_unsensitivity(drug_sen_df, useful_cell, 
                            p=get_sen_unsen_p, c=get_sen_unsen_c, indicator='auc')
    diff_cell = get_sen_unsen_res['sensitivity_cell_ls']+get_sen_unsen_res['unsensitivity_cell_ls']
    diff_rpkm_expr_mat = rpkm_expr_mat[get_sen_unsen_res['sensitivity_cell_ls']+get_sen_unsen_res['unsensitivity_cell_ls']]

    # wgcna
    print("""

    =============
    ||  wgcna  ||
    =============
    
    """)
    
    WGCNA_analysis_res = WGCNA_analysis(rpkm_expr_mat, drug_sen_df, n_top_mad=2000, plot_flag=True)
    # WGCNA_signi_ls = list(WGCNA_analysis_res['cluster_df'].loc['turquoise']['gene'])

    # limma
    print("""

    =============
    ||  limma  ||
    =============

    """)

    condition_table = np.array(['sen']*len(get_sen_unsen_res['sensitivity_cell_ls'])+
                        ['unsen']*len(get_sen_unsen_res['unsensitivity_cell_ls']))
    do_limma_res = do_limma(expr_mat=diff_rpkm_expr_mat, 
            condition_table=condition_table, 
            logFC_th=1, pvalue_th=0.05)
    # get_diff_from_limma_res = get_diff_from_limma(do_limma_res['limma_output'], logFC_th=2, pvalue_th=1e-7)


    # cor
    print("""

    ============
    ||  corr  ||
    ============
    
    """)

    cor_res_df = cor_onebyone(expr_mat=rpkm_expr_mat, 
                drug_sen_array=np.array(drug_sen_df['auc']))
    # significent_gene_ls = get_significient_from_cor(cor_res_df, r2_th=0.25, pvalue_th=1e-14)['significent_gene_ls']



    # interaction

    # len(set(WGCNA_signi_ls) & set(significent_gene_ls) & set(get_diff_from_limma_res['diff_gene_ls']))

    return {
        'useful_cell':useful_cell,
        'get_sen_unsen_res':get_sen_unsen_res,
        'WGCNA_analysis_res':WGCNA_analysis_res,
        'do_limma_res':do_limma_res,
        'cor_res_df':cor_res_df
    }


def pre_post_processing_diff_corr_wgcna():
    # pre
    #drug_sen_path = '/y/Bondi/data/VIZOME/nature_aml_drug_sen.tsv'
    #all_drug_sen_df = pd.read_csv(drug_sen_path, sep='\t', index_col=0)
    ch_drug_sen_path = '/y/Bondi/data/changhai_data/PDC体外用药数据汇总：持续更新ing20200722(1).xlsx'
    ch_drug_sen_df = load_changhai_drug_sen(ch_drug_sen_path)    

    rpkm_path = '/y/Bondi/data/nature_aml_fpkm.txt'
    rpkm_expr_mat = pd.read_csv(rpkm_path, sep='\t', index_col=0)

    ## there can do some filter to sample
    Variants_path =  '/y/Bondi/data/changhai_data/LAML_nature_PDC.maf'
    Variants_df = pd.read_csv(Variants_path, sep='\t')
    ### 获取 FLT3-ITD 的病人
    FLT3ITD_Variants_df = Variants_df.loc[Variants_df['Hugo_Symbol']=='FLT3'].loc[Variants_df['ITDorNOT']=='ITD']
    FLT3ITD_sample_ls = list(set(FLT3ITD_Variants_df['Tumor_Sample_Barcode']))
    FLT3ITD_useful_sample_ls = list(set(rpkm_expr_mat.columns) & set(FLT3ITD_sample_ls))


    # do
    do_diff_corr_wgcna_res = do_diff_corr_wgcna(drug_ls[drug_id], 
                    all_drug_sen_df=all_drug_sen_df,  ## drug sen data can be filter in the func
                    rpkm_expr_mat=rpkm_expr_mat[FLT3ITD_useful_sample_ls],  # there do the filter, sample is control by expr_mat
                    get_sen_unsen_p=0.2 , 
                    get_sen_unsen_c = None)
    
    
    # post
    ## wgcna
    do_diff_corr_wgcna_res['WGCNA_analysis_res']['moduleEigengenes_drug_cor_df'].sort_values('Adjusted R-squared')
    
    wgcna_ls_0 = do_diff_corr_wgcna_res['WGCNA_analysis_res']['cluster_df'].loc['green']['gene']
    
    print(len(wgcna_ls_0))

    for g in wgcna_ls_0:
        print(g)

    ## limma
    do_diff_corr_wgcna_res['do_limma_res']['limma_output'].sort_values('adj.P.Val')
    
    diff_limma_gene_ls_0 = get_diff_from_limma(
        do_diff_corr_wgcna_res['do_limma_res']['limma_output'], logFC_th=1.5, pvalue_th=1)['diff_gene_ls']

    print(len(diff_limma_gene_ls_0))

    for g in diff_limma_gene_ls_0:
        print(g)

    ## corr
    do_diff_corr_wgcna_res['cor_res_df'].sort_values('adj_r_squared')

    sign_cor_gene_ls_0 = get_significient_from_cor(do_diff_corr_wgcna_res['cor_res_df'], r2_th=0.1, pvalue_th=0.05)['significent_gene_ls']

    print(len(sign_cor_gene_ls_0))

    for g in sign_cor_gene_ls_0:
        print(g)


#####################
# mutation analysis #
#####################
# -------- mut -------
# VIZOME
def bar_plot_fusion_drug_sen(Variants_path, drug_sen_path, drug, plot_gene_ls, text, specific_sample_ls=None):
    '''
    input:
        Variants_path
        drug_sen_path
        drug
        text
        specific_sample_ls
    '''
    drug_sen_df = pd.read_csv(drug_sen_path, sep='\t', index_col=0)
    Variants_df = pd.read_csv(Variants_path, sep='\t')

    drug_sen_sample_ls = list(set(drug_sen_df.loc[drug]['lab_id']))


    ##
    mut_cell_sen_dict = {'gene':[], 'cell line':[], 'drug sensitivity':[]}
    for gene in plot_gene_ls:
        ## get gene_mut cell list
        specific_mut_sample_ls = Variants_df.loc[Variants_df['symbol']==gene]['labId'] # <-
        specific_mut_sample_ls = list(set(specific_mut_sample_ls))

        if specific_sample_ls is not None:
            useful_sample_ls = list(set(specific_mut_sample_ls)&set(drug_sen_sample_ls)&set(specific_sample_ls))
        else:
            useful_sample_ls = list(set(specific_mut_sample_ls)&set(drug_sen_sample_ls))
        # print(len(specific_sample_ls))
        ## get specific cell line drug sen
        # for cell_line in specific_sample_ls:
        mut_cell_sen_dict['drug sensitivity'] += list(drug_sen_df.loc[drug].set_index(['lab_id']).loc[useful_sample_ls]['auc']) # <-
        mut_cell_sen_dict['gene'] += [gene]*len(useful_sample_ls)
        mut_cell_sen_dict['cell line'] += useful_sample_ls

    mut_cell_sen_df = pd.DataFrame(mut_cell_sen_dict).set_index('gene')

    '''
    top_genes : a list of genes with top number of mut
    mut_cell_sen_df : 一个dataframe，描述每种mut有哪些cell line，这些cell line在我们的药下ic50 or AUC是多少
    '''
    avg_drug_sen_ls = []
    gene_ls_for_plot = []
    ls_for_plot = []
    for gene in plot_gene_ls:
        try:
            # print(list(mut_cell_sen_df.loc[[gene]]['drug sensitivity']))
            ls_for_plot.append(list(mut_cell_sen_df.loc[[gene]]['drug sensitivity']))
            avg_drug_sen_ls.append(np.mean(list(mut_cell_sen_df.loc[[gene]]['drug sensitivity'])))
            gene_ls_for_plot.append(gene)
        except:
            continue

    # ===        
    def sort_by_list(sort_by_ls, to_sort_ls):
        sort_index = np.array(sort_by_ls).argsort()
        sorted_ls = [to_sort_ls[i] for i in sort_index]
        return sorted_ls

    gene_ls_for_plot = sort_by_list(sort_by_ls=avg_drug_sen_ls, to_sort_ls=gene_ls_for_plot)
    ls_for_plot = sort_by_list(sort_by_ls=avg_drug_sen_ls, to_sort_ls=ls_for_plot) # !! avg_drug_sen_ls只拿来排序，作图用不到

    x = np.arange(len(gene_ls_for_plot))
    y = list(map(np.mean, ls_for_plot))
    std_err = [np.std(a, ddof = 1) for a in ls_for_plot]

    n = len(gene_ls_for_plot)
    x = np.arange(len(gene_ls_for_plot))
    y = list(map(np.mean, ls_for_plot))
    std_err = [np.std(a, ddof = 1)/math.sqrt(n) for a in ls_for_plot]


    scatter_x = []
    scatter_y = []
    for i_d in range(0,len(ls_for_plot)): ###################
        for c in ls_for_plot[i_d]:
            #x = np.random.normal(i, 0.02)
            scatter_x.append(np.random.normal(i_d, 0.1))
            #for_plot_x.append(i_d)
            scatter_y.append(c)





    import matplotlib.pyplot as plt
    '''
    x=[1,2,3,4,5]
    #数据集
    y=[20,44,21,64,46]
    #误差列表
    std_err=[1,2,5,3,2]
    '''
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
    ax.set_xlabel('common mutation ( with FLT3ITD mut )',font2)
    ax.set_ylabel('AUC',font2)
    
    ax.scatter(scatter_x, scatter_y,marker='*',color='r',s=1, zorder=1)

    plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
            rotation_mode="anchor")
    #显示图形
    ax.set_title('%s'%drug)
    plt.show()

## 还要选出BCOR CEBPA top back 的sample
def bar_plot_fusion_drug_sen(Variants_path, drug_sen_path, drug, plot_gene_ls, text, specific_sample_ls=None):
    '''
    input:
        Variants_path
        drug_sen_path
        drug
        text
        specific_sample_ls
    '''
    drug_sen_df = pd.read_csv(drug_sen_path, sep='\t', index_col=0)
    Variants_df = pd.read_csv(Variants_path, sep='\t')

    drug_sen_sample_ls = list(set(drug_sen_df.loc[drug]['lab_id']))


    ##
    mut_cell_sen_dict = {'gene':[], 'cell line':[], 'drug sensitivity':[]}
    for gene in plot_gene_ls:
        ## get gene_mut cell list
        specific_mut_sample_ls = Variants_df.loc[Variants_df['symbol']==gene]['labId'] # <-
        specific_mut_sample_ls = list(set(specific_mut_sample_ls))

        if specific_sample_ls is not None:
            useful_sample_ls = list(set(specific_mut_sample_ls)&set(drug_sen_sample_ls)&set(specific_sample_ls))
        else:
            useful_sample_ls = list(set(specific_mut_sample_ls)&set(drug_sen_sample_ls))
        # print(len(specific_sample_ls))
        ## get specific cell line drug sen
        # for cell_line in specific_sample_ls:
        mut_cell_sen_dict['drug sensitivity'] += list(drug_sen_df.loc[drug].set_index(['lab_id']).loc[useful_sample_ls]['auc']) # <-
        mut_cell_sen_dict['gene'] += [gene]*len(useful_sample_ls)
        mut_cell_sen_dict['cell line'] += useful_sample_ls

    mut_cell_sen_df = pd.DataFrame(mut_cell_sen_dict).set_index('gene')

    '''
    top_genes : a list of genes with top number of mut
    mut_cell_sen_df : 一个dataframe，描述每种mut有哪些cell line，这些cell line在我们的药下ic50 or AUC是多少
    '''
    avg_drug_sen_ls = []
    gene_ls_for_plot = []
    ls_for_plot = []
    for gene in plot_gene_ls:
        try:
            # print(list(mut_cell_sen_df.loc[[gene]]['drug sensitivity']))
            ls_for_plot.append(list(mut_cell_sen_df.loc[[gene]]['drug sensitivity']))
            avg_drug_sen_ls.append(np.mean(list(mut_cell_sen_df.loc[[gene]]['drug sensitivity'])))
            gene_ls_for_plot.append(gene)
        except:
            continue

    # ===        
    def sort_by_list(sort_by_ls, to_sort_ls):
        sort_index = np.array(sort_by_ls).argsort()
        sorted_ls = [to_sort_ls[i] for i in sort_index]
        return sorted_ls

    gene_ls_for_plot = sort_by_list(sort_by_ls=avg_drug_sen_ls, to_sort_ls=gene_ls_for_plot)
    #
    top3 = gene_ls_for_plot[-3:]
    back3 = gene_ls_for_plot[:3]
    if 'BCOR' in back3 or 'CEBPA' in top3:
        #

        ls_for_plot = sort_by_list(sort_by_ls=avg_drug_sen_ls, to_sort_ls=ls_for_plot) # !! avg_drug_sen_ls只拿来排序，作图用不到

        x = np.arange(len(gene_ls_for_plot))
        y = list(map(np.mean, ls_for_plot))
        std_err = [np.std(a, ddof = 1) for a in ls_for_plot]

        n = len(gene_ls_for_plot)
        x = np.arange(len(gene_ls_for_plot))
        y = list(map(np.mean, ls_for_plot))
        std_err = [np.std(a, ddof = 1)/math.sqrt(n) for a in ls_for_plot]


        scatter_x = []
        scatter_y = []
        for i_d in range(0,len(ls_for_plot)): ###################
            for c in ls_for_plot[i_d]:
                #x = np.random.normal(i, 0.02)
                scatter_x.append(np.random.normal(i_d, 0.1))
                #for_plot_x.append(i_d)
                scatter_y.append(c)





        import matplotlib.pyplot as plt
        '''
        x=[1,2,3,4,5]
        #数据集
        y=[20,44,21,64,46]
        #误差列表
        std_err=[1,2,5,3,2]
        '''
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
        
        ax.scatter(scatter_x, scatter_y,marker='*',color='r',s=1, zorder=1)

        plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
                rotation_mode="anchor")
        #显示图形
        ax.set_title('%s'%drug)
        plt.show()

        return 1
    return 0


# Changhai
def bar_plot_fusion_drug_sen(Variants_path, drug_sen_path, drug, plot_gene_ls, text, specific_sample_ls=None):
    '''
    input:
        Variants_path
        drug_sen_path
        drug
        text
        specific_sample_ls
    '''
    drug_sen_df = pd.read_csv(drug_sen_path, sep='\t', index_col=0)
    Variants_df = pd.read_csv(Variants_path, sep='\t')

    drug_sen_sample_ls = list(set(drug_sen_df.loc[drug]['lab_id']))


    ##
    mut_cell_sen_dict = {'gene':[], 'cell line':[], 'drug sensitivity':[]}
    for gene in plot_gene_ls:
        ## get gene_mut cell list
        specific_mut_sample_ls = Variants_df.loc[Variants_df['Hugo_Symbol']==gene]['Tumor_Sample_Barcode']
        specific_mut_sample_ls = list(set(specific_mut_sample_ls))

        if specific_sample_ls is not None:
            useful_sample_ls = list(set(specific_mut_sample_ls)&set(drug_sen_sample_ls)&set(specific_sample_ls))
        else:
            useful_sample_ls = list(set(specific_mut_sample_ls)&set(drug_sen_sample_ls))
        # print(len(specific_sample_ls))
        ## get specific cell line drug sen
        # for cell_line in specific_sample_ls:
        mut_cell_sen_dict['drug sensitivity'] += list(drug_sen_df.loc[drug].set_index(['lab_id']).loc[useful_sample_ls]['auc'])
        mut_cell_sen_dict['gene'] += [gene]*len(useful_sample_ls)
        mut_cell_sen_dict['cell line'] += useful_sample_ls

    mut_cell_sen_df = pd.DataFrame(mut_cell_sen_dict).set_index('gene')

    '''
    top_genes : a list of genes with top number of mut
    mut_cell_sen_df : 一个dataframe，描述每种mut有哪些cell line，这些cell line在我们的药下ic50 or AUC是多少
    '''
    avg_drug_sen_ls = []
    gene_ls_for_plot = []
    ls_for_plot = []
    for gene in plot_gene_ls:
        try:
            # print(list(mut_cell_sen_df.loc[[gene]]['drug sensitivity']))
            ls_for_plot.append(list(mut_cell_sen_df.loc[[gene]]['drug sensitivity']))
            avg_drug_sen_ls.append(np.mean(list(mut_cell_sen_df.loc[[gene]]['drug sensitivity'])))
            gene_ls_for_plot.append(gene)
        except:
            continue

    # ===        
    def sort_by_list(sort_by_ls, to_sort_ls):
        sort_index = np.array(sort_by_ls).argsort()
        sorted_ls = [to_sort_ls[i] for i in sort_index]
        return sorted_ls

    gene_ls_for_plot = sort_by_list(sort_by_ls=avg_drug_sen_ls, to_sort_ls=gene_ls_for_plot)
    ls_for_plot = sort_by_list(sort_by_ls=avg_drug_sen_ls, to_sort_ls=ls_for_plot) # !! avg_drug_sen_ls只拿来排序，作图用不到

    x = np.arange(len(gene_ls_for_plot))
    y = list(map(np.mean, ls_for_plot))
    std_err = [np.std(a, ddof = 1) for a in ls_for_plot]

    n = len(gene_ls_for_plot)
    x = np.arange(len(gene_ls_for_plot))
    y = list(map(np.mean, ls_for_plot))
    std_err = [np.std(a, ddof = 1)/math.sqrt(n) for a in ls_for_plot]


    scatter_x = []
    scatter_y = []
    for i_d in range(0,len(ls_for_plot)): ###################
        for c in ls_for_plot[i_d]:
            #x = np.random.normal(i, 0.02)
            scatter_x.append(np.random.normal(i_d, 0.1))
            #for_plot_x.append(i_d)
            scatter_y.append(c)





    import matplotlib.pyplot as plt
    '''
    x=[1,2,3,4,5]
    #数据集
    y=[20,44,21,64,46]
    #误差列表
    std_err=[1,2,5,3,2]
    '''
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
    ax.set_xlabel('common mutation ( with FLT3ITD mut )',font2)
    ax.set_ylabel('AUC',font2)
    
    ax.scatter(scatter_x, scatter_y,marker='*',color='r',s=1, zorder=1)

    plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
            rotation_mode="anchor")
    #显示图形
    ax.set_title('%s'%drug)
    plt.show()





# -------- Fusion -------


def bar_plot_fusion_drug_sen(clinical_annotation_path, drug_sen_path, drug, text, specific_sample_ls=None):
    '''
    import:
        clinical_annotation_path
        drug_sen_path
        drug
    '''
    # clinical_annotation_path = '/y/Bondi/data/changhai_data/nature_aml_clinical_annotation.txt'
    clinical_annotation_df = pd.read_csv(clinical_annotation_path, sep='\t', encoding='gb18030', index_col=0)

    # drug_sen_path = '/y/Bondi/data/VIZOME/nature_aml_drug_sen.tsv'
    all_drug_sen_df = pd.read_csv(drug_sen_path, sep='\t', index_col=0)

    ##
    # drug = 'Venetoclax'
    ##

    useful_cell = list( set(clinical_annotation_df['finalFusion'].index) & set(all_drug_sen_df.loc[drug]['lab_id']) )
    if specific_sample_ls is not None:
        useful_cell = list(set(useful_cell)&set(specific_sample_ls))
    
    useful_drug_sen_df = all_drug_sen_df.loc[drug].set_index(['lab_id']).loc[useful_cell]
    useful_clinical_annotation_df = clinical_annotation_df.loc[useful_cell]


    # 找出各个fusion，有哪些cell line， 找出这些cellline的 drug sen
    sample_Series = pd.Series(useful_drug_sen_df.index)
    sample_Series.index = useful_drug_sen_df.index
    
    fusion_cell_sen_df = pd.concat(
        [ useful_clinical_annotation_df[['finalFusion']], useful_drug_sen_df, sample_Series], axis=1
        ).set_index('finalFusion')

    # prepare for plot
    fusion_ls_for_plot = list(set(fusion_cell_sen_df.index))
    drug_sen_ls_for_plot = []
    avg_drug_sen_ls = []
    for fusion in fusion_ls_for_plot:
        drug_sen_ls_for_plot.append( list(fusion_cell_sen_df.loc[[fusion]]['auc']) )
        avg_drug_sen_ls.append( np.mean(list(fusion_cell_sen_df.loc[[fusion]]['auc'])) )

    ## 从大到小
    fusion_ls_for_plot = sort_by_list(sort_by_ls=avg_drug_sen_ls, to_sort_ls=fusion_ls_for_plot)
    drug_sen_ls_for_plot = sort_by_list(sort_by_ls=avg_drug_sen_ls, to_sort_ls=drug_sen_ls_for_plot)
    # avg_drug_sen_ls = sorted(avg_drug_sen_ls) # avg_drug_sen_ls只拿来排序，作图用不到


    n = len(fusion_ls_for_plot)
    x = np.arange(len(fusion_ls_for_plot))
    y = list(map(np.mean, drug_sen_ls_for_plot))
    '''
    !! 小心：如果某种fusion只有一个样本的数据，这里会warning，并且算不出std
    '''
    std_err = [np.std(a, ddof = 1)/math.sqrt(n) for a in drug_sen_ls_for_plot]

    scatter_x = []
    scatter_y = []
    for i_d in range(0,len(drug_sen_ls_for_plot)): ###################
        for c in drug_sen_ls_for_plot[i_d]:
            #x = np.random.normal(i, 0.02)
            scatter_x.append(np.random.normal(i_d, 0.1))
            #for_plot_x.append(i_d)
            scatter_y.append(c)

    # plot
    error_params=dict(elinewidth=2,ecolor='black',capsize=3)#设置误差标记参数
    font2 = {
            'weight' : 'normal',
            'size' : 15,
            }

    fig, ax = plt.subplots(figsize=[5,3], facecolor='white')

    ax.bar(x,y,yerr=std_err,error_kw=error_params, zorder=0)
    ax.set_xticks(np.arange(len(fusion_ls_for_plot)))
    ax.set_xticklabels(fusion_ls_for_plot)
    ax.set_xlabel('Final Fusion',font2)
    ax.set_ylabel('AUC',font2)
        
    ax.scatter(scatter_x, scatter_y,marker='*',color='r',s=1, zorder=1)

    plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
            rotation_mode="anchor")
    ax.set_title('%s %s'%(text,drug))

    plt.show()

    return {
        'useful_cell':useful_cell,
        'fusion_cell_sen_df':fusion_cell_sen_df,
        'fusion_ls_for_plot':fusion_ls_for_plot,
        'avg_drug_sen_ls':avg_drug_sen_ls
    }


####################
# subtype analysis #
####################




#############
# load data #
#############

def load_changhai_drug_sen(ch_drug_sen_path):
    # ch_drug_sen_path = '/y/Bondi/data/changhai_data/PDC体外用药数据汇总：持续更新ing20200722(1).xlsx'
    ch_drug_sen_df = pd.read_excel(ch_drug_sen_path, header=2, index_col=0).iloc[:,3:]
    return ch_drug_sen_df



####################
# boxplot expr~mut #
####################


def for_test(expr_gene, mut_gene, Variants_df, rpkm_expr_mat):
    #drug = 'Venetoclax'
    #expr_gene = 'MYCN'
    #mut_gene = 'CEBPA'

    # pre
    #Variants_path =  '/y/Bondi/data/changhai_data/LAML_nature_PDC.maf'
    #Variants_df = pd.read_csv(Variants_path, sep='\t')

    #drug_sen_path = '/y/Bondi/data/VIZOME/nature_aml_drug_sen.tsv'
    #all_drug_sen_df = pd.read_csv(drug_sen_path, sep='\t', index_col=0)

    #rpkm_path = '/y/Bondi/data/nature_aml_fpkm.txt'
    #rpkm_expr_mat = pd.read_csv(rpkm_path, sep='\t', index_col=0)



    # do
    #drug_sen_df = all_drug_sen_df.loc[drug].set_index(['lab_id']) # can change

    #useful_cell = set(rpkm_expr_mat.columns) & set(drug_sen_df.index) # can change


    #rpkm_expr_mat = rpkm_expr_mat[useful_cell]
    rpkm_expr_mat = clean_data(rpkm_expr_mat, dup='avg')

    #drug_sen_df = drug_sen_df.loc[useful_cell]

    ## specific mut
    Specifical_Mut_Variants_df = Variants_df.loc[Variants_df['Hugo_Symbol']==mut_gene]
    Specifical_Mut_sample_ls = list(set(Specifical_Mut_Variants_df['Tumor_Sample_Barcode']))

    Specifical_Mut_useful_sample_ls = list(set(rpkm_expr_mat.columns) & set(Specifical_Mut_sample_ls))
    specifical_Wild_useful_sample_ls = list(rpkm_expr_mat.columns.difference(Specifical_Mut_useful_sample_ls))


    Specifical_mut_expr_Series = rpkm_expr_mat.loc[expr_gene, Specifical_Mut_useful_sample_ls]
    Specifical_wild_expr_Series = rpkm_expr_mat.loc[expr_gene, specifical_Wild_useful_sample_ls]


    ## for plot
    scatter_y = np.array(list(Specifical_mut_expr_Series) + list(Specifical_wild_expr_Series))
    scatter_x = np.array(len(Specifical_mut_expr_Series)*[1]+len(Specifical_wild_expr_Series)*[2])
    rnd_for_add = np.array([np.random.normal(0, 0.05) for i in range(len(scatter_x))])
    scatter_x = scatter_x + rnd_for_add

    ## plot
    font2 = {
            'weight' : 'normal',
            'size' : 15,
            }

    fig, ax = plt.subplots()

    ax.boxplot(x=[Specifical_mut_expr_Series, Specifical_wild_expr_Series], zorder=1)
    ax.scatter(scatter_x, scatter_y, marker='*',color='g',s=1, zorder=0)
    ax.set_xticks(np.array([1,2]))
    ax.set_xticklabels(['Mut', 'WT'], font2)
    ax.set_xlabel('Mut %s'%mut_gene, font2)
    ax.set_ylabel('Expression(RPKM) %s'%expr_gene, font2)

    plt.show()


    %R -i Specifical_mut_expr_Series
    %R -i Specifical_wild_expr_Series
    %R t_res = t.test(Specifical_mut_expr_Series,Specifical_wild_expr_Series,paired = FALSE) # 这里就先省略方差齐性检验
    %R print(t_res)


def for_test2(expr_gene1, expr_gene2, rpkm_expr_mat, test = None):

    gene1_expr_Series = rpkm_expr_mat.loc[expr_gene1]
    gene2_expr_Series = rpkm_expr_mat.loc[expr_gene2]
    

    do_cor(x=gene1_expr_Series, 
            y=gene2_expr_Series, 
            x_name='Expression(RPKM) %s'%expr_gene1, 
            y_name='Expression(RPKM) %s'%expr_gene2, 
            c=None, 
            plot=True, 
            figsize=None, 
            plot_text=test)
    '''
    # plot
    font2 = {
            'weight' : 'normal',
            'size' : 13,
            }

    fig, ax = plt.subplots()

    plt.scatter(gene1_expr_Series, gene2_expr_Series)
    ax.set_xlabel('Expression(RPKM) %s'%expr_gene1, font2)
    ax.set_ylabel('Expression(RPKM) %s'%expr_gene2, font2)

    plt.show()

    %R -i gene1_expr_Series
    %R -i gene1_expr_Series
    %R 
    '''

