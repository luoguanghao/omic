library(tidyverse)
library(pheatmap)
library(ConsensusClusterPlus)

path_perk = c('EIF2AK3','EIF2S1','ATF4')
path_ire1 = c('ERN1','MAP3K5','MAP2K7','MAPK8')
path_atf6 = c('ATF6','HSPA5')

urp_list = c(pat_perk=pat_perk,path_ire1=path_ire1,path_atf6=path_atf6)

exp_mat = read.csv('~/my_project/data/VIZOME/nature_aml_log2_fpkm_clean.txt',sep='\t')


annotation_col = data.frame(
  Cluster = factor(c(rep('perk',length(path_perk)) , rep('ire1',length(path_ire1)) , rep('atf6',length(path_atf6)) ))
)
rownames(annotation_col) = c(path_perk, path_ire1, path_atf6)

pheatmap(
    exp_mat[c(path_perk,path_ire1,path_atf6),],
    cluster_cols = TRUE,
    cluster_rows = FALSE,
    scale = 'row',
    color = colorRampPalette(c("green","black", "red"))(50),
    filename='./upr_hmap.pdf',
    cellwidth = 2, cellheight = 15,
    annotation_row = annotation_col,
    border_color = "white"
    )


mat_for_cc = exp_mat[c(path_perk,path_ire1,path_atf6),]
results = ConsensusClusterPlus(t(scale(t(mat_for_cc))),maxK=8)


### 这里是取4类的结果来分析
class_df = as.data.frame(results[[4]][['consensusClass']],col.names=c('class'))
colnames(class_df)=c('class')
class_df['sample'] = rownames(class_df)

tmp1 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['1']]
tmp2 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['2']]
tmp3 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['3']]
tmp4 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['4']]

#### 4c
# annotation_col = data.frame(
#  Cluster_col = factor(c(rep('c1',length(tmp1)) , rep('c2',length(tmp2)) , rep('c3',length(tmp3)) , rep('c4',length(tmp4))))
# )
# rownames(annotation_col) = c(tmp1,tmp2,tmp3,tmp4)

annotation_row = data.frame(
  Cluster_row = factor(c(rep('perk',length(path_perk)) , rep('ire1',length(path_ire1)) , rep('atf6',length(path_atf6)) ))
)
rownames(annotation_row) = c(path_perk, path_ire1, path_atf6)

anno_colors<-list(
  Cluster_col=c(c1="red",c2="blue",c3="green",c4="pink"),
  Cluster_row=c(perk="red",ire1="blue",atf6="green")
)


pheatmap(
    exp_mat[c(path_perk,path_ire1,path_atf6),],
    cluster_cols = TRUE,
    cluster_rows = FALSE,
    scale = 'row',
    color = colorRampPalette(c("green","black", "red"))(50),
    filename='./upr_hmap.pdf',
    cellwidth = 2, cellheight = 15,
    annotation_row = annotation_row,
    annotation_col = annotation_col,
    border_color = "white",
    #color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    annotation_colors = anno_colors
    )










