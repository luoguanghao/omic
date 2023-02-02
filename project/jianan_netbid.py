import pandas as pd
import numpy as np
import math, csv, pickle
import random, os
import matplotlib.pyplot as plt

%load_ext rpy2.ipython
%R library(NetBID2)

def make_name_R_useful(name_array, symb_be_change='-',symb_change_to='.', front_add = 'X'):
    new_name_array = np.array([front_add+item.replace(symb_be_change,symb_change_to) for item in name_array])
    return new_name_array

def get_sensitivity_unsensitivity(GDSC_df, cell_list, p=None, c=None, indicator='IC50'):
    '''
    取得药物敏感性排头和拍尾的细胞系名称
    对于药物敏感性，是IC50小的敏感性强， IC50大的敏感性不强
    GDSC_df: 
    cell IC50 AUC
    -------------
    xx   8.8  88
    
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



def DO_NetBID_pipeline():

    %R -i exp_mat
    %R cal_eset <- generate.eset(exp_mat = exp_mat, phenotype_info = NULL,feature_info = NULL, annotation_info = "exp mat from /y/Bondi/other_method/NetBID/AML/data/nodup_nature_aml_fpkm_SJARACNe.tsv")

    %R analysis.par = list()
    %R analysis.par$cal.eset = cal_eset




    # phenotype prepare
    G1 = make_name_R_useful(np.array(get_sensitivity_unsensitivity_res['sensitivity_cell_ls']))
    G0 = make_name_R_useful(np.array(get_sensitivity_unsensitivity_res['unsensitivity_cell_ls']))
    %R -i G1
    %R -i G0

    # DE
    # %R DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='sensitivity',G0_name='unsensitivity')
    %R DE_gene_bid_limma <- getDE.limma.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='sensitivity',G0_name='unsensitivity')

    # get network
    tf_network_path = '/y/Bondi/other_method/NetBID/AML/SJAR/output_tf_sjaracne_project_2021-2-20_out_.final/consensus_network_ncol_.txt'
    sig_network_path = '/y/Bondi/other_method/NetBID/AML/SJAR/output_sig_sjaracne_project_2021-2-20_out_.final/consensus_network_ncol_.txt'

    %R -i tf_network_path
    %R -i sig_network_path
    %R analysis.par$tf.network.file = tf_network_path
    %R analysis.par$sig.network.file = sig_network_path
    %R analysis.par$tf.network <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
    %R analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)

    # Merge network first
    %R analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)

    # Get activity matrix 
    # !!!!!
    %R ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')

    # generate eset

    %R analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=NULL,feature_info=NULL,annotation_info='activity in net-dataset tf&sig')

    # DA
    # %R DE_gene_bid_limma <- getDE.limma.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='sensitivity',G0_name='unsensitivity')
    %R DA_driver_bid_limma <- getDE.limma.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='sensitivity',G0_name='unsensitivity')

    # post DA
    %R analysis.par$DE <- list()
    %R analysis.par$DA <- list()
    %R comp_name <- 'sensitivity.Vs.unsensitivity' # Each comparison must has a name

    # Save comparison result to list element in analysis.par, with comparison name
    %R analysis.par$DE[[comp_name]] <- DE_gene_bid_limma
    %R analysis.par$DA[[comp_name]] <- DA_driver_bid_limma

    %R DE_gene_comb <- combineDE(DE_list=list('sensitivity.Vs.unsensitivity'=analysis.par$DE$`sensitivity.Vs.unsensitivity`,'sensitivity.Vs.unsensitivity'=analysis.par$DE$`sensitivity.Vs.unsensitivity`))
    %R DA_driver_comb <- combineDE(DE_list=list('sensitivity.Vs.unsensitivity'=analysis.par$DA$`sensitivity.Vs.unsensitivity`,'sensitivity.Vs.unsensitivity'=analysis.par$DA$`sensitivity.Vs.unsensitivity`))

    ## PLOT !!{
    %R draw.combineDE(DE_gene_comb)

    %R draw.combineDE(DA_driver_comb)
    ## }

    ## Check DataFrame {
    DA_driver_bid_limma = %R DA_driver_bid_limma
    DE_gene_bid_limma = %R DE_gene_bid_limma

    DA_driver_bid_limma.sort_values('logFC').iloc[30:40]
    DA_driver_bid_limma.sort_values('P.Value')
    DE_gene_bid_limma.sort_values('logFC', ascending=False)
    ## }

    # save master table
    %R all_comp <- names(analysis.par$DE)
    %R db.preload(use_level='gene',use_spe='human',update=FALSE)

    %R analysis.par$final_ms_tab <- generate.masterTable(use_comp=all_comp, DE=analysis.par$DE, DA=analysis.par$DA, target_list=analysis.par$merge.network$target_list, tf_sigs=tf_sigs,z_col='Z-statistics', display_col=c('logFC','P.Value'), main_id_type='external_gene_name')

    # Adv. analysis
    %R ms_tab <- analysis.par$final_ms_tab ## get the master table data frame
    %R ms_tab <- ms_tab[which(ms_tab$Size>=30 & ms_tab$Size <=1000),] 
    %R comp_name <- 'sensitivity.Vs.unsensitivity' ## get the comparison name


    ## volcanoPlot
    %R sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',logFC_col=sprintf('logFC.%s_DA',comp_name),Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.1,Pv_thre=1e-8,main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,label_cex = 1)

    %R sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',logFC_col=sprintf('logFC.%s_DA',comp_name),Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.1,Pv_thre=1e-2,main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,label_cex = 1)
    ## 


    # Get the DE data frame of target genes
    %R DE <- analysis.par$DE[[comp_name]]
    %R driver_list <- rownames(sig_driver) # The rownames is the originalID_label

    ## GSEA plot {
    %R draw.GSEA.NetBID(DE=DE,profile_col='logFC',profile_trend='neg2pos',name_col='ID',driver_list = driver_list,show_label=ms_tab[driver_list,'gene_label'],driver_DA_Z=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],driver_DE_Z=ms_tab[driver_list,sprintf('Z.%s_DE',comp_name)],target_list=analysis.par$merge.network$target_list,top_driver_number=60,target_nrow=2,target_col='RdBu', left_annotation = 'high in others',right_annotation = 'high in G4',main=comp_name,target_col_type='DE',Z_sig_thre=1.64,profile_sig_thre = 1.64)
    ## }



def do():
    # get data    
    exp_mat_path = '/y/Bondi/other_method/NetBID/AML/data/nodup_nature_aml_fpkm_SJARACNe.tsv'
    exp_mat = pd.read_csv(exp_mat_path, sep='\t', index_col=0).drop('geneSymbol',axis=1)

    drug_sen_path = '/y/Bondi/data/VIZOME/nature_aml_drug_sen.tsv'
    all_drug_sen_df = pd.read_csv(drug_sen_path, sep='\t', index_col=0)

    # drug_ls = ['GDC-0941', 'BEZ235', 'PI-103', 'Venetoclax', 'Sorafenib']
    drug = 'Venetoclax'

    drug_sen_df = all_drug_sen_df.loc['Venetoclax'].set_index('lab_id')
    useful_sample = set(drug_sen_df.index) & set(exp_mat.columns)
    drug_sen_df = drug_sen_df.loc[useful_sample]
    exp_mat = exp_mat[useful_sample]

    get_sensitivity_unsensitivity_res = get_sensitivity_unsensitivity(drug_sen_df, drug_sen_df.index, p=0.2, indicator='ic50')


































