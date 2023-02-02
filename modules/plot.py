# only for jupyter
%load_ext rpy2.ipython
%R library(maftools)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def get_data():
    pass


def sort_by_list(sort_by_ls, to_sort_ls):
    sort_index = np.array(sort_by_ls).argsort()
    sorted_ls = [to_sort_ls[i] for i in sort_index]
    return sorted_ls

def plot_jianan_graph():
    ##########
    # get data
    get_data_res = get_data()

    gene_maf_df = get_data_res['gene_maf_df']
    drug_sen_df = get_data_res['changhai_drug_sen']
    maf_obj = get_data_res['maf_obj']

    #################
    # get top N genes
    top_n = 10
    #
    GeneSummary_df = %R getGeneSummary(maf_obj)
    top_genes = GeneSummary_df.iloc[:top_n]['Hugo_Symbol']

    ##########################################################
    # 找出各个突变，有哪些cell line， 找出这些cellline的 drug sen
    '''
    changhai_gene_maf_df
    changhai_drug_sen_df
    top_genes
    '''
    drug = 'BEZ235 (NVP-BEZ235, Dactolisib)'
    #
    mut_cell_sen_dict = {'gene':[], 'cell line':[], 'drug sensitivity':[]}

    for gene in top_genes:
        for cell_line in gene_maf_df.loc[gene]['Tumor_Sample_Barcode']:
            try:
                if drug_sen_df[drug][cell_line] == drug_sen_df[drug][cell_line]:
                    mut_cell_sen_dict['gene'].append(gene)
                    mut_cell_sen_dict['cell line'].append(cell_line)
                    mut_cell_sen_dict['drug sensitivity'].append(drug_sen_df[drug][cell_line])
            except:
                print('!!! %s'%cell_line)

    mut_cell_sen_df = pd.DataFrame(mut_cell_sen_dict).set_index('gene')


    ##############
    # plot boxplot
    '''
    top_genes : a list of genes with top number of mut
    mut_cell_sen_df : 一个dataframe，描述每种mut有哪些cell line，这些cell line在我们的药下ic50 or AUC是多少
    '''
    # prepare for plot
    ls_for_plot = []
    avg_drug_sen_ls = []
    gene_ls_for_plot = []
    for gene in top_genes:
        try:
            # print(list(mut_cell_sen_df.loc[[gene]]['drug sensitivity']))
            ls_for_plot.append(list(mut_cell_sen_df.loc[[gene]]['drug sensitivity']))
            avg_drug_sen_ls.append(np.mean(list(mut_cell_sen_df.loc[[gene]]['drug sensitivity'])))
            gene_ls_for_plot.append(gene)
        except:
            continue
    ## 从大到小
    gene_ls_for_plot = sort_by_list(sort_by_ls=avg_drug_sen_ls, to_sort_ls=gene_ls_for_plot)
    ls_for_plot = sort_by_list(sort_by_ls=avg_drug_sen_ls, to_sort_ls=ls_for_plot)

    # plot
    fig, ax = plt.subplots()

    ax.boxplot(ls_for_plot)
    ax.set_xticks(np.arange(len(gene_ls_for_plot)))
    ax.set_xticklabels(gene_ls_for_plot)
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
            rotation_mode="anchor")
    plt.show()






#####













































































