options(repr.matrix.max.cols=50, repr.matrix.max.rows=300) 

library(tidyverse)
library(stringr)

# ...
a = "FAM63A
FAM63B
FAM188A
FAM188B
ALG13
ATXN3
ATXN3L
FAM105A
JOSD1
JOSD2
OTUB1
OTUB2
OTUD1
OTUD3
OTUD4
OTUD5
OTUD6A
OTUD6B
OTUD7A
OTUD7B
OTULIN
TNFAIP3
VCPIP1
YOD1
ZRANB1
BAP1
UCHL1
UCHL3
UCHL5
CYLD
PAN2
USP1
USP2
USP3
USP4
USP5
USP6
USP7
USP8
USP9X
USP9Y
USP10
USP11
USP12
USP13
USP14
USP15
USP16
USP17L2
USP18
USP19
USP20
USP21
USP22
USP24
USP25
USP26
USP27X
USP28
USP29
USP30
USP31
USP32
USP33
USP38
ZUFSP
"
de_ubiquitin_ls = strsplit(a,'\n')[[1]]

# LAML
exp_dir = '~/my_project/data/TCGA/LAML/TCGA-LAML.htseq_fpkm.tsv.gz'
probeMap_dir = '~/my_project/data/TCGA/LAML/gencode.v22.annotation.gene.probeMap'

#exp_dir = '~/my_project/data/TCGA/DLBC/TCGA-DLBC.htseq_fpkm.tsv.gz'
#probeMap_dir = '~/my_project/data/TCGA/DLBC/gencode.v22.annotation.gene.probeMap'

##import TCGA data
exp_mat = read.csv(exp_dir, sep='\t')
probeMap_df = read.csv(probeMap_dir, sep='\t')

probeMap_df <- probeMap_df %>%  dplyr::select(1,2)
colnames(probeMap_df) <- c("Ensembl_ID","gene_name")

exp_mat <- inner_join(probeMap_df, exp_mat, by = "Ensembl_ID") %>% dplyr::select(-1)
exp_mat <- aggregate(.~ gene_name, exp_mat, mean)
exp_mat <- exp_mat %>% column_to_rownames("gene_name")




#####################################
# 所有别的基因和query gene进行cor计算 #
# res_df_1
#####################################
# input ,delete the gene with 0 var
t_exp_mat = as.data.frame(t(exp_mat))
t_exp_mat = t_exp_mat[,apply(t_exp_mat,2,var)!=0] # 删去方差为0的基因

query_gene = 'USP7' ######
cor_threshold = 0.1 #####


#cor_res = cor(dplyr::select(t_exp_mat, -c(all_of(query_gene))), t_exp_mat[,query_gene])

#cor_res = na.omit( as.data.frame(cor_res) )
#colnames(cor_res) = 'correlation'
#cor_res$gene = rownames(cor_res)
#cor_res = cor_res[,c('gene','correlation')]
#cor_res = cor_res[order(cor_res$correlation),]

#cor_cor_res = cor_res%>%filter(abs(correlation)>cor_threshold)

df = dplyr::select(t_exp_mat, -c(all_of(query_gene)))
x = t_exp_mat[,query_gene]
y_i_ls = colnames(df)

res_df = list(gene=c(),r=c(),p=c())
##
for(y_i in y_i_ls){
    y = df[[y_i]]
    rrr = cor.test( x=x, y=y )
    res_df$gene = c(res_df$gene, y_i)
    res_df$r = c(res_df$r, rrr$estimate[[1]])
    res_df$p = c(res_df$p, rrr$p.value)
}
res_df = as.data.frame(res_df)

# na.omit(res_df[order(res_df$r),])%>%filter(abs(r)>cor_threshold)

res_df_1 = na.omit(res_df[order(res_df$r),])

#########################################
# 筛选出阈值基因，然后和de_niqui基因做corr #
# res_df_2
#########################################

sel_genes_1 = (res_df_1%>%filter(r>cor_threshold))$gene
# sel_genes_1 = res_df_1[order(res_df_1$r, decreasing=TRUE),]$gene[1:3000]

res_df_2 = cor(t_exp_mat[,sel_genes_1], t_exp_mat[,de_ubiquitin_ls])

res_df_2 = as.data.frame(res_df_2)

##########################################################
# 筛选出基因:其与其他de_uniqu基因的相关性全都不大于usp7的基因 #
# res_df_2[sel_rows,]
##########################################################


sel_rows = c()
for(i_r in 1:dim(res_df_2)[1]){
    if(sum(res_df_2[i_r,]>res_df_2$USP7[i_r])==0){
        sel_rows = c(sel_rows,i_r)
    }
}
res_df_2[sel_rows,]


## 格式的整理
res_df_3 = res_df_2[sel_rows,][,colnames(res_df_2[sel_rows,])!='USP7']
res_df_3$USP7 = res_df_2[sel_rows,]$USP7
res_df_3[  ,c(colnames(res_df_3)[dim(res_df_3)[2]], colnames(res_df_3)[-dim(res_df_3)[2]])  ]



write.csv(res_df_3[  ,c(colnames(res_df_3)[dim(res_df_3)[2]], colnames(res_df_3)[-dim(res_df_3)[2]])  ],
          file='LAML_sel_gene_matrix_th0.1.csv',
          row.names=TRUE)




#=====================================================
# OLD
#=====================================================
# cor() 的结果是相关系数
#
#


library(tidyverse)
library(stringr)


query_gene = 'USP7'
cor_threshold = 0.5



exp_dir = '/y/Archive/Bondi/data/TCGA/LAML/TCGA-LAML.htseq_fpkm.tsv.gz'
probeMap_dir = '/y/Archive/Bondi/data/TCGA/LAML/gencode.v22.annotation.gene.probeMap'

##import TCGA data
exp_mat = read.csv(exp_dir, sep='\t')
probeMap_df = read.csv(probeMap_dir, sep='\t')

probeMap_df <- probeMap_df %>%  dplyr::select(1,2)
colnames(probeMap_df) <- c("Ensembl_ID","gene_name")

exp_mat <- inner_join(probeMap_df, exp_mat, by = "Ensembl_ID") %>% dplyr::select(-1)
exp_mat <- aggregate(.~ gene_name, exp_mat, mean)
exp_mat <- exp_mat %>% column_to_rownames("gene_name")


t_exp_mat = as.data.frame(t(exp_mat))

cor_res = cor(dplyr::select(t_exp_mat, -c(query_gene)), t_exp_mat[,query_gene]) # 有一些warning!!
# Note: Using an external vector in selections is ambiguous.
# ℹ Use `all_of(query_gene)` instead of `query_gene` to silence this message.
# ℹ See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
# This message is displayed once per session.

Warning message in cor(dplyr::select(t_exp_mat, -c(query_gene)), t_exp_mat[, query_gene]):
“the standard deviation is zero”

cor_res = na.omit( as.data.frame(cor_res) )
colnames(cor_res) = 'correlation'
cor_res$gene = rownames(cor_res)
cor_res = cor_res[,c('gene','correlation')]
cor_res = cor_res[order(cor_res$correlation),]

cor_cor_res = cor_res%>%filter(abs(correlation)>cor_threshold)


## for循环计算 ##

df = dplyr::select(t_exp_mat, -c(all_of(query_gene)))
x = t_exp_mat[,query_gene]
y_i_ls = colnames(df)

res_df = list(gene=c(),r=c(),p=c())
##
for(y_i in y_i_ls){
    y = df[[y_i]]
    rrr = cor.test( x=x, y=y )
    res_df$gene = c(res_df$gene, y_i)
    res_df$r = c(res_df$r, rrr$estimate[[1]])
    res_df$p = c(res_df$p, rrr$p.value)
}
res_df = as.data.frame(res_df)

################













if (0){

for (i in 1:nrow(data)){
  for (r in i:nrow(data)){
    g1=rownames(data)[i]
    g2=rownames(data)[r]
    c_r=cor(as.numeric(data[i,]),as.numeric(data[r,]),method="pearson")
    p=cor.test(as.numeric(data[i,]),as.numeric(data[r,]),method ="pearson")[[3]]
    ##保存每一步的数据，而不可直接以空向量作为每一步运行的结果
    gene_name1=c(gene_name1,g1)
    gene_name2=c(gene_name2,g2)
    cor_r=c(cor_r,c_r)
    pvalue=c(pvalue,p)
       }
}


}






















































