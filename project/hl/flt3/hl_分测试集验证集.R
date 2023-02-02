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


sp_drug_sen_df = drug_sen_df%>%filter(inhibitor==drug)
## 看看药敏的分布
ggplot(data=drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],]) + 
  geom_histogram(aes(x=auc),bins=30)

df_dot_plot = drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],][order(drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],]$auc),]
df_dot_plot$lab_id = factor(df_dot_plot$lab_id,level=as.vector(df_dot_plot$lab_id))

plot(x=1:423,y=df_dot_plot$ic50)
ggplot(data=df_dot_plot) + geom_point(aes(x=lab_id, y=auc)) + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))




drug = our_drugs[1]
# analysis
sp_drug_sen_df = drug_sen_df%>%filter(inhibitor==drug)
## divide
### get sen unsen
top_bottom_number = round(dim(sp_drug_sen_df)[1]*0.3)
sp_drug_sen_df = drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],][order(drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],]$auc),]
sp_drug_sen_df$lab_id = factor(sp_drug_sen_df$lab_id,level=as.vector(sp_drug_sen_df$lab_id))
rownames(sp_drug_sen_df) = sp_drug_sen_df$lab_id
sp_sample = intersect(sp_drug_sen_df$lab_id,colnames(exp_mat))

sen_sample = (sp_drug_sen_df[sp_drug_sen_df$lab_id %in% sp_sample,] %>% filter(auc<182.3626))$lab_id ### should change
unsen_sample = (sp_drug_sen_df[sp_drug_sen_df$lab_id %in% sp_sample,] %>% filter(auc>232.7781))$lab_id ### should change

### divide
nn = 0.8 ### change!!
sub<-sample(1:length(sen_sample),round(length(sen_sample)*nn))
t_sen_sample = as.vector(sen_sample[sub])
v_sen_sample = as.vector(sen_sample[-sub])

sub<-sample(1:length(unsen_sample),round(length(unsen_sample)*nn))
t_unsen_sample = as.vector(unsen_sample[sub])
v_unsen_sample = as.vector(unsen_sample[-sub])


## select feature

condition_table = c(rep('sen',length(t_sen_sample)),rep('unsen',length(t_unsen_sample)))
exp_mat_for_anal = exp_mat[,c(t_sen_sample,t_unsen_sample)]
useful_gene = rownames(as.data.frame(apply(exp_mat_for_anal,1,var)) %>% filter(apply(exp_mat_for_anal, 1, var)>0))
exp_mat_for_anal = exp_mat_for_anal[useful_gene,]


### do limma
design <- model.matrix(~condition_table)
colnames(design) <- levels(condition_table)
rownames(design) <- colnames(exp_mat_for_anal)
fit <- lmFit(exp_mat_for_anal, design)
fit <- eBayes(fit, trend=TRUE)
output <- topTable(fit, coef=2,n=Inf)
output = output[order(output$P.Value),]

## post analysis
### hl feature
ND_gene_ls_dir = '/y/Archive/Bondi/jupyter/hl/3-16-hl-根据ND基因聚类两类/hl-ND-gene.txt'
ND_gene_ls <- read.table(ND_gene_ls_dir,header = FALSE,sep = '\n')


intersect(rownames(output%>%filter(adj.P.Val<0.05)), as.vector(ND_gene_ls[['V1']])) # look intersect
intersect(rownames(output%>%filter(P.Value<0.05)), as.vector(ND_gene_ls[['V1']])) # look intersect
exp_mat[as.vector(ND_gene_ls[['V1']]),] # look intersect
output[as.vector(ND_gene_ls[['V1']]),] # look intersect

### pheatmap
anno_nr_df = data.frame(type=c(rep('sen',length(t_sen_sample)),rep('unsen',length(t_unsen_sample))),row.names=c(t_sen_sample,t_unsen_sample))
pheatmap(
    exp_mat_for_anal[as.vector(ND_gene_ls[['V1']]),],
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    scale = 'row',
    color = colorRampPalette(c("green","black", "red"))(50),
    filename='./cencluster_hmap1.png',
    cellwidth = 8, cellheight = 18,border=FALSE,
    cutree_cols=2, annotation_col  = anno_nr_df
    )
plot(x=1:134, y=sp_drug_sen_df[colnames(exp_mat_for_anal),]$auc)
plot(x=1:134, y=sp_drug_sen_df[colnames(exp_mat_for_anal),]$ic50)

