import pandas as pd
import numpy as np
import math, csv, pickle
import random, os

%load_ext rpy2.ipython


def prepare_for_deseq2_from_TCGA(expr_mat):
    condition_table = []
    for i in expr_mat.columns:
        condition_table.append(distinguish_cancer(i))
    condition_table = np.array(condition_table)

    colData =  pd.DataFrame({'Condition':condition_table},index=list(expr_mat.columns))

    return {"condition_table":condition_table, 'colData':colData}


def prepare_for_deseq2():
    pass
    #return {"condition_table":condition_table, 'colData':colData}

def do_deseq2(expr_mat, condition_table, colData):
    '''
    expr_mat : genes*samples
    condition_table：1,0的array
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



def differential_expression_from_TCGA(expr_mat, dump=False, plot_for_check=False, up_down='all'):
    '''
    up_down: up/down/all
    '''
    ans_prepare_for_deseq2 = prepare_for_deseq2(expr_mat)
    condition_table = ans_prepare_for_deseq2['condition_table']
    colData = ans_prepare_for_deseq2['colData']

    ans_do_deseq2 = do_deseq2(expr_mat, condition_table, colData)
    DESeq_res = ans_do_deseq2['DESeq_res']
    expr_norm = pd.DataFrame(ans_do_deseq2['expr_norm'], columns=expr_count.columns, index=expr_count.index)
    
    if dump:
        with open("DESeq2_result.pkl","wb") as file:
            pickle.dump(ans_do_deseq2, file, True)
    
    gene_ls = np.array(expr_mat.index)
    
    if up_down == 'up':
        diff_genes = list(set(gene_ls[DESeq_res['log2FoldChange']>2]) & set(gene_ls[DESeq_res['pvalue'] < 0.05]))
    elif up_down == 'down':
        diff_genes = list(set(gene_ls[DESeq_res['log2FoldChange']<-2]) & set(gene_ls[DESeq_res['pvalue'] < 0.05]))
    else:
        diff_genes = list(set(gene_ls[abs(DESeq_res['log2FoldChange'])>2]) & set(gene_ls[DESeq_res['pvalue'] < 0.05]))

    return {'diff_genes':diff_genes, 'DESeq_res':DESeq_res,'expr_norm':expr_norm}


def CPM2Count(expr_mat):
    pass




















































