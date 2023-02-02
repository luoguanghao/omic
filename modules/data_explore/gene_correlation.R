options(repr.matrix.max.cols=50, repr.matrix.max.rows=300) 

library(tidyverse)
library(stringr)



#####################################
# 所有别的基因和query gene进行cor计算 #
# res_df_1
#####################################
# query_gene 想拿来做相关的基因
# exp_mat 表达矩阵，也可能是GSVA的矩阵：genes*samples
# 
gene_cor_with_other <- function(query_gene, exp_mat){
    t_exp_mat = as.data.frame(t(exp_mat))
    t_exp_mat = t_exp_mat[,apply(t_exp_mat,2,var)!=0] # 删去方差为0的基因

    # query_gene = 'USP7' ######
    # cor_threshold = 0.1 #####

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
    return(res_df_1)
}

# 等会归位
# 选择一些与usp7相关性强的基因，
# 然后再将与usp7相关性强的基因与去泛素化基因做相关
# 筛选一些usp7相关基因：其与去泛素化基因
# 
#
usp7_pipeline <- function(query_gene,de_ubiquitin_ls,exp_mat,filename=NULL){

    res_df_1 <- gene_cor_with_other(query_gene, exp_mat)

    #########################################
    # 筛选出阈值基因，然后和de_niqui基因做corr #
    # res_df_2
    #########################################

    sel_genes_1 = (res_df_1%>%filter(r>cor_threshold))$gene
    # sel_genes_1 = res_df_1[order(res_df_1$r, decreasing=TRUE),]$gene[1:3000]

    res_df_2 = cor(t_exp_mat[,sel_genes_1], t_exp_mat[,de_ubiquitin_ls])

    res_df_2 = as.data.frame(res_df_2) # row是usp7相关基因，col是去泛素化基因

    ##########################################################
    # 筛选出基因:其与其他de_uniqu基因的相关性全都不大于usp7的基因 #
    # res_df_2[sel_rows,]
    ##########################################################

    sel_rows = c()
    for(i_r in 1:dim(res_df_2)[1]){
        if(sum(res_df_2[i_r,]>res_df_2$USP7[i_r])==0){ # 比较某一个基因与其它去泛素化基因的相关性是否不大于usp7
            sel_rows = c(sel_rows,i_r)
        }
    }
    # res_df_2[sel_rows,]

    ## 格式的整理
    res_df_3 = res_df_2[sel_rows,][,colnames(res_df_2[sel_rows,])!='USP7'] # 不要usp7那一列
    res_df_3$USP7 = res_df_2[sel_rows,]$USP7 # usp7放最后
    res_df_3[  ,c(colnames(res_df_3)[dim(res_df_3)[2]], colnames(res_df_3)[-dim(res_df_3)[2]])  ]# usp7放最前

    if(is.null(filename)==FALSE){
        write.csv(res_df_3[  ,c(colnames(res_df_3)[dim(res_df_3)[2]], colnames(res_df_3)[-dim(res_df_3)[2]])  ],
                file=filename, row.names=TRUE)
    }
    return(list(res_df_1=res_df_1, res_df_2=res_df_2, res_df_3=res_df_3))
}























