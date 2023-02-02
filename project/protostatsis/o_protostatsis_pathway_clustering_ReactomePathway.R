library(tidyverse)
library(DESeq2)
library(GSVA)
library(pheatmap)
library(ConsensusClusterPlus)

deseq2_normalization <- function(deseq2_exp_mat, condition_table){

    # =========
    dds <- DESeqDataSetFromMatrix(deseq2_exp_mat, 
                                DataFrame(condition_table), 
                                design= ~ condition_table)
    dds <- dds[rowSums(counts(dds)) > 1,]
    dds2 <- DESeq(dds)
    resultsNames(dds2)
    # acquire the results using function results(), and assign to res
    res <- results(dds2)
    # view the summary of results
    summary(res)
    # resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
    # output
    #head(resdata)

    deseq2_normalized_exp_mat = as.data.frame(counts(dds2,normalize=TRUE))

    return(deseq2_normalized_exp_mat)
}


pw_list_path = './data/pw_list.txt'
pw_list_df = read_tsv(pw_list_path)

# ==== TCGA data preprocess
exp_path = '~/my_project/data/TCGA/BRCA/TCGA-BRCA.htseq_counts.tsv.gz'
probeMap_df = read.csv('~/my_project/data/TCGA/BRCA/gencode.v22.annotation.gene.probeMap',sep='\t')
###
probeMap_df <- probeMap_df %>%  dplyr::select(1,2)
colnames(probeMap_df) <- c("Ensembl_ID","gene_name")

exp_mat <- inner_join(probeMap_df, exp_mat, by = "Ensembl_ID") %>% dplyr::select(-1)
exp_mat <- aggregate(.~ gene_name, exp_mat, mean)
exp_mat <- exp_mat %>% column_to_rownames("gene_name")

# =========== TCGA deseq2 normalization


deseq2_exp_mat = exp_mat
condition_table = factor(c(rep('sen',608),rep('unsen',609)))
# deseq2_exp_mat = as.matrix(exp_mat%>%column_to_rownames('gene_id'))
deseq2_exp_mat = apply(deseq2_exp_mat,2,as.integer)
rownames(deseq2_exp_mat) = rownames(exp_mat)

deseq2_normalized_exp_mat = deseq2_normalization(deseq2_exp_mat,condition_table)




# === select pathway gene
exp_mat_cc = deseq2_normalized_exp_mat[unique(pw_list_df$F),] # 选择pathway基因
exp_mat_cc = na.omit(exp_mat_cc)

# Consensus Cluster #################
# input of ConsensusClusterPlus: gene×sample ,this func is do the cluster for sample
##
results = ConsensusClusterPlus(t(scale(t(exp_mat_cc))),maxK=8)
save(results, file = "./result/consenCluster_res_vizome.RData")


## post cluster analysis

# pw_sort_ls = read.table('./pathway_sort_ls.txt',header = FALSE,sep = '\n')

### 这里是取6类的结果来分析
class_df = as.data.frame(results[[6]][['consensusClass']],col.names=c('class'))
colnames(class_df)=c('class')
class_df['sample'] = rownames(class_df)

tmp1 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['1']]
tmp2 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['2']]
tmp3 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['3']]
tmp4 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['4']]
tmp5 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['5']]
tmp6 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['6']]

annotation_col = data.frame(
  Cluster = factor(c(rep('c1',length(tmp1)) , rep('c2',length(tmp2)) , rep('c3',length(tmp3)) , rep('c4',length(tmp4)), rep('c5',length(tmp5)), rep('c6',length(tmp6))))
)
rownames(annotation_col) = c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)


gene_ls = rownames(exp_mat_cc)

anno_df = list(gene=c(), anno=c())
for(g in gene_ls){
    anno_df$gene = c(anno_df$gene, g)
    anno_df$anno = c(anno_df$anno, paste(unique((pw_list_df%>%filter(F==g))$A), collapse='_'))
}
anno_df = as.data.frame(anno_df)%>%filter(gene%in%rownames(exp_mat_cc))
anno_df = na.omit(anno_df[order(anno_df$anno),])


annotation_row$Cluster_row[annotation_row$Cluster_row=='']='_'

annotation_col = data.frame(
  Cluster_col = factor(c(rep('c1',length(tmp1)) , rep('c2',length(tmp2)) , rep('c3',length(tmp3)) , rep('c4',length(tmp4)), rep('c5',length(tmp5)), rep('c6',length(tmp6))))
)
rownames(annotation_col) = c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)

annotation_row = data.frame(
    Cluster_row = as.vector(anno_df$anno)
)
annotation_row$Cluster_row[annotation_row$Cluster_row=='']='N'
annotation_row$Cluster_row = factor(annotation_row$Cluster_row)
rownames(annotation_row) = anno_df$gene


anno_colors<-list(
  Cluster_col=c(c1="red",c2="blue",c3="green",c4="pink",c5='black',c6='yellow'),
  Cluster_row=c(CMA="red",UPR="blue",UPS="green",UPS_CMA="pink",N='black')
)


pheatmap(
    (exp_mat_cc[anno_df$gene,c(tmp1,tmp2,tmp3,tmp4)]),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    scale = 'row',
    color = colorRampPalette(c("blue","white", "red"))(50),
    filename='./cencluster_hmap3.pdf',
    cellwidth = 2, cellheight = 10,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = anno_colors[1],
    border_color = "white"
    )























































































