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



# get data
## mut data

maf_df = read.csv('/y/Archive/Bondi/data/VIZOME/LAML_nature_PDC.maf',sep='\t')
ITD_sample_ls = as.vector((maf_df %>% filter(Hugo_Symbol=='FLT3'&ITDorNOT=='ITD'))['Tumor_Sample_Barcode'][['Tumor_Sample_Barcode']])
ITD_sample_ls = str_replace_all(paste('X',ITD_sample_ls,sep=''),'-','.')



## exp data
exp_mat = read.csv('/y/Archive/Bondi/data/VIZOME/nature_aml_log2_fpkm.txt',sep='\t')
# 处理重复基因，滤出需要的基因
exp_mat = aggregate(.~ Symbol, exp_mat, mean)
exp_mat = exp_mat %>% column_to_rownames('Symbol')
sample_ls = colnames(exp_mat)


## drug sen data
drug_sen_df = read.csv('/y/Archive/Bondi/data/VIZOME/nature_aml_drug_sen.tsv',sep='\t')
drug_sen_df$lab_id = lapply(drug_sen_df$lab_id, function(x) str_replace_all(paste('X',x,sep=''),'-','.'))
drug_sen_df$lab_id = as.factor(as.character(drug_sen_df$lab_id))
drug_ls = drug_sen_df[['inhibitor']][!duplicated(drug_sen_df[['inhibitor']])]

our_drugs = c('Midostaurin','Quizartinib (AC220)','Crenolanib','Gilteritinib (ASP-2215)')

# get sp drug 
cat(our_drugs[1]) ## change
sp_drug_sen_df = drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],][order(drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],]$auc),]
sp_drug_sen_df$lab_id = factor(sp_drug_sen_df$lab_id,level=as.vector(sp_drug_sen_df$lab_id))
rownames(sp_drug_sen_df) = sp_drug_sen_df$lab_id

sp_sample = intersect(intersect(sp_drug_sen_df$lab_id,colnames(exp_mat)),ITD_sample_ls)

# get sen unsen
top_bottom_number = round(dim(sp_drug_sen_df %>% filter(lab_id %in% ITD_sample_ls))[1]*0.3) # set the sen unsen ####
sen_cutoff = (sp_drug_sen_df %>% filter(lab_id %in% ITD_sample_ls))[top_bottom_number,]$auc
unsen_cutoff = (sp_drug_sen_df %>% filter(lab_id %in% ITD_sample_ls))[
    dim((sp_drug_sen_df %>% filter(lab_id %in% ITD_sample_ls)))[1]-top_bottom_number,]$auc

sen_sample = as.vector((sp_drug_sen_df[sp_drug_sen_df$lab_id %in% sp_sample,] %>% filter(auc<sen_cutoff))$lab_id )
unsen_sample = as.vector((sp_drug_sen_df[sp_drug_sen_df$lab_id %in% sp_sample,] %>% filter(auc>unsen_cutoff))$lab_id )

# divide
nn = 0.7  #### set the discovery coho ####

sub<-sample(1:length(sen_sample),round(length(sen_sample)*nn))
t_sen_sample = as.vector(sen_sample[sub])
v_sen_sample = as.vector(sen_sample[-sub])

sub<-sample(1:length(unsen_sample),round(length(unsen_sample)*nn))
t_unsen_sample = as.vector(unsen_sample[sub])
v_unsen_sample = as.vector(unsen_sample[-sub])

# set pheno
phenotype_info = as.data.frame(list(sample=c(t_sen_sample,t_unsen_sample), 
                      sen_unsen=c(rep('sen',length(t_sen_sample)),rep('unsen',length(t_unsen_sample)))))
rownames(phenotype_info) = phenotype_info$sample


# NetBID #########################################
# construct analysis.par 
exp_mat_nb = exp_mat
G0 = t_sen_sample
G1 = t_unsen_sample
G0_name = 'sensitivity'
G1_name = 'unsensitivity'
comp_name <- 'sensitivity.Vs.unsensitivity' # Each comparison must has a name

tf_network_path = '/y/Archive/Bondi/other_method/NetBID/AML/SJAR/output_tf_sjaracne_project_2021-2-20_out_.final/consensus_network_ncol_.txt'
sig_network_path = '/y/Archive/Bondi/other_method/NetBID/AML/SJAR/output_sig_sjaracne_project_2021-2-20_out_.final/consensus_network_ncol_.txt'
##
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
















