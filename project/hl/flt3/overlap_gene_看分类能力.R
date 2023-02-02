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



# input data

## mut data
maf_df = read.csv('~/my_project/data/VIZOME/LAML_nature_PDC.maf',sep='\t')
ITD_sample_ls = as.vector((maf_df %>% filter(Hugo_Symbol=='FLT3'&ITDorNOT=='ITD'))['Tumor_Sample_Barcode'][['Tumor_Sample_Barcode']])
ITD_sample_ls = str_replace_all(paste('X',ITD_sample_ls,sep=''),'-','.')

## exp data
exp_mat = read.csv('~/my_project/data/VIZOME/nature_aml_log2_fpkm.txt',sep='\t')
# 处理重复基因，滤出需要的基因
exp_mat = aggregate(.~ Symbol, exp_mat, mean)
exp_mat = exp_mat %>% column_to_rownames('Symbol')
sample_ls = colnames(exp_mat)

## drug sen data
drug_sen_df = read.csv('~/my_project/data/VIZOME/nature_aml_drug_sen.tsv',sep='\t')
drug_sen_df$lab_id = lapply(drug_sen_df$lab_id, function(x) str_replace_all(paste('X',x,sep=''),'-','.'))
drug_sen_df$lab_id = as.factor(as.character(drug_sen_df$lab_id))
drug_ls = drug_sen_df[['inhibitor']][!duplicated(drug_sen_df[['inhibitor']])]


our_drugs = c('Midostaurin','Quizartinib (AC220)','Crenolanib','Gilteritinib (ASP-2215)')




# ND gene input and do overlap

ND_gene_dir = './NDgene.txt'
ND_gene_file = read_tsv(ND_gene_dir)

MIDO_ls = na.omit(ND_gene_file$MIDO)
QUIZ_ls = na.omit(ND_gene_file$QUIZ)
GILTER_ls = na.omit(ND_gene_file$GILTER)
CRENO_ls = na.omit(ND_gene_file$CRENO)

all_gene = unique(c(MIDO_ls,QUIZ_ls,GILTER_ls,CRENO_ls))
length(all_gene)

df_for_plot = as.data.frame(list( gene=rep(all_gene,4),
                                drug=c(rep('MIDO',length(all_gene)),rep('QUIZ',length(all_gene)),rep('GILTER',length(all_gene)),rep('CRENO',length(all_gene))),
                                value=c(all_gene%in%MIDO_ls,all_gene%in%QUIZ_ls,all_gene%in%GILTER_ls,all_gene%in%CRENO_ls)
                                ))



ND_df_mean=group_by(df_for_plot%>%filter(value==TRUE), gene) %>% summarise(count = n())
df_for_plot$gene = factor(df_for_plot$gene, levels = as.vector(ND_df_mean[order(ND_df_mean$count),][['gene']]))
ND_overlap_gene = as.vector( (ND_df_mean%>%filter(count==4))$gene )

cols=c('TRUE'='red','FALSE'='black')

df_for_plot%>%ggplot(aes(x=drug,y=gene))+
  geom_tile(aes(fill=value),color="white",size=1)+ #color和size分别指定方块边线的颜色和粗细
  scale_x_discrete("",expand = c(0,0))+ #不显示横纵轴的label文本；画板不延长
  scale_y_discrete("",expand = c(0,0))+
  coord_equal(0.1) +
  scale_fill_manual(values = cols)+ #指定自定义的颜色
  theme(
    axis.text.x.bottom = element_text(size=10),axis.text.y.left = element_text(size = 12), #修改坐标轴文本大小
    axis.ticks = element_blank(), #不显示坐标轴刻度
    legend.title = element_blank() #不显示图例title
  )
ggsave("tmp.pdf",device = "pdf",width = 9,height = 25)


# 看看overlap基因的热图

# ND_overlap_gene = as.vector( (df_mean%>%filter(vars==4))$gene )

### pheatmap
# anno_nr_df = data.frame(type=c(rep('sen',length(t_sen_sample)),rep('unsen',length(t_unsen_sample))),row.names=c(t_sen_sample,t_unsen_sample))
pheatmap(
    exp_mat[ND_overlap_gene,itd_exp_drug_sample],
    #cluster_cols = FALSE,
    #cluster_rows = FALSE,
    scale = 'row',
    color = colorRampPalette(c("green","black", "red"))(50),
    #filename='./cencluster_hmap1.png',
    cellwidth = 8, cellheight = 18,border=FALSE,
    #cutree_cols=2, 
    #annotation_col  = anno_nr_df
    )


# sp drug drug 1



## 看看药敏的分布
ggplot(data=drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],]) + 
  geom_histogram(aes(x=auc),bins=30)

ggplot(data=drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],]) + 
  geom_histogram(aes(x=ic50),bins=30)

df_dot_plot = drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],][order(drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],]$auc),]
df_dot_plot$lab_id = factor(df_dot_plot$lab_id,level=as.vector(df_dot_plot$lab_id))

## 看AUC排名与ic50的关系
plot(x=1:423,y=df_dot_plot$ic50)
ggplot(data=df_dot_plot) + geom_point(aes(x=lab_id, y=auc)) + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))



## choose sp sample
itd_exp_drug_sample = intersect(intersect(ITD_sample_ls,colnames(exp_mat)), (drug_sen_df%>%filter(inhibitor==our_drugs[1]))$lab_id)

## divide
### get sen unsen

sp_drug_sen_df = drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],]%>%filter(lab_id%in%itd_exp_drug_sample)
sp_drug_sen_df = sp_drug_sen_df[order(sp_drug_sen_df$auc),]
sp_drug_sen_df$lab_id = factor(sp_drug_sen_df$lab_id,level=as.vector(sp_drug_sen_df$lab_id))

top_bottom_number = round(dim(sp_drug_sen_df)[1]*0.3)

sen_sample = as.vector(sp_drug_sen_df[1:top_bottom_number,]$lab_id)
unsen_sample = as.vector(sp_drug_sen_df[(dim(sp_drug_sen_df)[1]-top_bottom_number):dim(sp_drug_sen_df)[1],]$lab_id)

# sp drug drug 1



## 看看药敏的分布
ggplot(data=drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],]) + 
  geom_histogram(aes(x=auc),bins=30)

ggplot(data=drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],]) + 
  geom_histogram(aes(x=ic50),bins=30)

df_dot_plot = drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],][order(drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],]$auc),]
df_dot_plot$lab_id = factor(df_dot_plot$lab_id,level=as.vector(df_dot_plot$lab_id))

## 看AUC排名与ic50的关系
plot(x=1:dim(df_dot_plot)[1],y=df_dot_plot$ic50)
ggplot(data=df_dot_plot) + geom_point(aes(x=lab_id, y=auc)) + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))



## choose sp sample
itd_exp_drug_sample = intersect(intersect(ITD_sample_ls,colnames(exp_mat)), (drug_sen_df%>%filter(inhibitor==our_drugs[1]))$lab_id)

## divide
### get sen unsen

sp_drug_sen_df = drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],]%>%filter(lab_id%in%itd_exp_drug_sample)
sp_drug_sen_df = sp_drug_sen_df[order(sp_drug_sen_df$auc),]
sp_drug_sen_df$lab_id = factor(sp_drug_sen_df$lab_id,level=as.vector(sp_drug_sen_df$lab_id))

top_bottom_number = round(dim(sp_drug_sen_df)[1]*0.3)

sen_sample = as.vector(sp_drug_sen_df[1:top_bottom_number,]$lab_id)
unsen_sample = as.vector(sp_drug_sen_df[(dim(sp_drug_sen_df)[1]-top_bottom_number):dim(sp_drug_sen_df)[1],]$lab_id)

## ssGSEA打分

ND_gene_ls = list(ND_overlap_gene=ND_overlap_gene)

gsva_matrix <- gsva(as.matrix(exp_mat[,c(sen_sample,unsen_sample)]), 
                    ND_gene_ls,
                    method='ssgsea', 
                    kcdf='Gaussian', 
                    abs.ranking=TRUE, 
                    parallel.sz=1)

df_plot = as.data.frame(t(gsva_matrix))

df_plot$sample = rownames(df_plot)

df_plot$sen_unsen = c(rep('sen',length(sen_sample)), rep('unsen',length(unsen_sample)))

ggplot(data=df_plot) + geom_boxplot(aes(x=sen_unsen, y=ND_overlap_gene,fill=sen_unsen)) + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))


































