library(stringr)
library(ggpubr) 
library(ggplot2)
library(GSEABase)
library(GSVAdata)
#data(c2BroadSets)
#c2BroadSets
 
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)
library(tidyverse)
library(pheatmap)
library(DESeq2)

library(edgeR)
library(limma)

library(NetBID2)



# NetBID #########################################
# example project\hl_分测试集验证集_netbid.py
# construct analysis.par 
exp_mat_nb = exp_mat # gene*sample
G0 = t_sen_sample # vector contain sample_name
G1 = t_unsen_sample
G0_name = 'sensitivity'
G1_name = 'unsensitivity'
comp_name <- 'sensitivity.Vs.unsensitivity' # Each comparison must has a name

tf_network_path = '/y/Archive/Bondi/other_method/NetBID/AML/SJAR/output_tf_sjaracne_project_2021-2-20_out_.final/consensus_network_ncol_.txt'
sig_network_path = '/y/Archive/Bondi/other_method/NetBID/AML/SJAR/output_sig_sjaracne_project_2021-2-20_out_.final/consensus_network_ncol_.txt'
#####

analysis.par = list()
analysis.par$cal.eset <- generate.eset(exp_mat = exp_mat_nb, phenotype_info = NULL,feature_info = NULL, annotation_info = "")


# DE
DE_gene_bid_limma <- getDE.limma.2G(eset=analysis.par$cal.eset,
                                    G0=G0, G1=G1, G0_name=G0_name, G1_name=G1_name)

# netbid

analysis.par$tf.network.file = tf_network_path
analysis.par$sig.network.file = sig_network_path
analysis.par$tf.network <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)

analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)

## cal activity
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')

analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,
                                            phenotype_info=NULL,feature_info=NULL,
                                            annotation_info='activity in net-dataset tf&sig')
DA_driver_bid_limma <- getDE.limma.2G(eset=analysis.par$merge.ac.eset,
                                      G0=t_sen_sample,G1=t_unsen_sample,G0_name='sensitivity',G1_name='unsensitivity')

## post analysis ##############################

analysis.par$DE <- list()
analysis.par$DA <- list()

#### comp_name <- 'sensitivity.Vs.unsensitivity' # Each comparison must has a name
analysis.par$DE[[comp_name]] <- DE_gene_bid_limma
analysis.par$DA[[comp_name]] <- DA_driver_bid_limma

DE_gene_comb <- combineDE( ## <<
    DE_list=list('sensitivity.Vs.unsensitivity'=analysis.par$DE$`sensitivity.Vs.unsensitivity`,'sensitivity.Vs.unsensitivity'=analysis.par$DE$`sensitivity.Vs.unsensitivity`))
DA_driver_comb <- combineDE(  # <<
    DE_list=list('sensitivity.Vs.unsensitivity'=analysis.par$DA$`sensitivity.Vs.unsensitivity`,'sensitivity.Vs.unsensitivity'=analysis.par$DA$`sensitivity.Vs.unsensitivity`))


all_comp <- names(analysis.par$DE)
db.preload(use_level='gene',use_spe='human',update=FALSE)
analysis.par$final_ms_tab <- generate.masterTable(use_comp=all_comp, DE=analysis.par$DE, 
                                                  DA=analysis.par$DA, target_list=analysis.par$merge.network$target_list, 
                                                  tf_sigs=tf_sigs,z_col='Z-statistics', display_col=c('logFC','P.Value'), 
                                                  main_id_type='external_gene_name')

ms_tab <- analysis.par$final_ms_tab ## get the master table data frame
ms_tab <- ms_tab[which(ms_tab$Size>=30 & ms_tab$Size <=1000),] 
#### comp_name <- 'sensitivity.Vs.unsensitivity' ## get the comparison name




## get driver
sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
                               logFC_col=sprintf('logFC.%s_DA',comp_name),Pv_col=sprintf('P.Value.%s_DA',comp_name),
                               logFC_thre=0.01,Pv_thre=1e-2,main=sprintf('Volcano Plot for %s_DA',comp_name),
                               show_label=TRUE,label_cex = 1)

## GSEA plot
DE <- analysis.par$DE[[comp_name]]
DA <- analysis.par$DA[[comp_name]]
driver_list <- rownames(sig_driver) # The rownames is the originalID_label
draw.GSEA.NetBID(DE=DE,profile_col='logFC',profile_trend='neg2pos',name_col='ID',
                 driver_list = driver_list,show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
                 driver_DE_Z=ms_tab[driver_list,sprintf('Z.%s_DE',comp_name)],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=60,target_nrow=2,target_col='RdBu', 
                 left_annotation = 'high in others',right_annotation = 'high in G4',main=comp_name,  ### <<<<
                 target_col_type='DE',Z_sig_thre=1.64,profile_sig_thre = 1.64)


## pheatmap
sig_driver_ls = sig_driver$gene_label

anno_nr_df = data.frame(type=c(rep('sen',length(t_sen_sample)),rep('unsen',length(t_unsen_sample))),row.names=c(t_sen_sample,t_unsen_sample))  # <<<<
pheatmap(ac_mat[sig_driver_ls,c(t_sen_sample,t_unsen_sample)],
        scale='row',
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        annotation_col  = anno_nr_df)
pheatmap(ac_mat[sig_driver_ls,c(t_sen_sample,t_unsen_sample)],
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        annotation_col  = anno_nr_df)
pheatmap(ac_mat[sig_driver_ls,c(t_sen_sample,t_unsen_sample)],
         scale='row',
        #cluster_cols = FALSE,
        #cluster_rows = FALSE,
        annotation_col  = anno_nr_df)

v_anno_nr_df = data.frame(type=c(rep('sen',length(v_sen_sample)),rep('unsen',length(v_unsen_sample))),row.names=c(v_sen_sample,v_unsen_sample))
pheatmap(ac_mat[sig_driver_ls,c(v_sen_sample,v_unsen_sample)],
         scale='row',
        #cluster_cols = FALSE,
        #cluster_rows = FALSE,
        annotation_col  = v_anno_nr_df
        )



















