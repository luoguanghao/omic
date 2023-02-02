library("clusterProfiler")
library(GSVA)

#############
# ssGSEA
#############
# do ssGSEA #
## input
# gsva_exp_mat:matrix ,rowname is gene, colnames is sample
# gene_ls: a list of gene
# other_information: df with 
#
## output
# a df with 1 col: ssGSEA score
do_ssGSEA <- function(gsva_exp_mat,gene_ls, other_information=NULL){
    gene_list = list(gene_ls=gene_ls)

    gsva_matrix <- gsva(gsva_exp_mat, 
                        gene_list,
                        method='ssgsea', 
                        kcdf='Gaussian', 
                        abs.ranking=TRUE, 
                        parallel.sz=1)

    gsva_res = data.frame(sample=colnames(gsva_matrix),ssGSEA_score=gsva_matrix[1,])

    return(gsva_res)
}


do_gsva <- function(exp_mat,pw_data,pw_g_list){
    # GSVA #################
    # exp_mat: gene×sample
    ##
    exp_mat = exp_mat[pw_data$Gene.symbol,] # 选择pathway基因
    exp_mat = na.omit(exp_mat)

    gsva_matrix <- gsva(as.matrix(exp_mat), 
                        pw_g_list,method='ssgsea', 
                        kcdf='Gaussian', 
                        abs.ranking=TRUE, parallel.sz=1)
    return(gsva_matrix)
}

# exp_mat : 表达矩阵：gene×sample
# gene_set : 基因集
# gene_set_type : 基因集的类型：支持gmt,fuscc,list,df等的格式
new_do_gsva <- function(exp_mat,gene_set,gene_set_type='gmt',min.sz=1,max.sz=Inf){
    # GSVA #################
    # exp_mat: gene×sample
    ##

    if(gene_set_type=='gmt'){
        pw_g_list = list()
        terms = unique(gene_set$term)
        for(t in terms){
            pw_g_list[[t]] = (gene_set%>%filter(term==t))$gene
        }
    }else if (gene_set_type=='fuscc') {
        gene_vec = c()
        pathway_vec = c()
        for(i_r in 1:dim(pw_data)[1]){
        # for(i_r in 1:3){
            gene = pw_data[i_r,1]
            pathways = as.character(unlist(strsplit(pw_data[i_r,2], split = "; ")))
            for(pw in pathways){
                gene_vec = c(gene_vec,gene)
                pathway_vec = c(pathway_vec,pw)
            }
        }
        gene_pathway_df = data.frame(gene=gene_vec, pathway=pathway_vec)
        gene_pathway_list <- split(as.matrix(gene_pathway_df)[,1], gene_pathway_df[,2])
        pw_g_list = split(as.matrix(gene_pathway_df)[,1], gene_pathway_df[,2])
    }else if (gene_set_type=='list') {
        pw_g_list=gene_set
    }else if (gene_set_type=='df') {
        pw_g_list = list()
        terms = unique(gene_set$term)
        for(t in terms){
            pw_g_list[[t]] = (gene_set%>%filter(term==t))$gene
        }
    }

    ##

    all_gene = c()
    for(i in 1:length(pw_g_list)){
        all_gene = c(all_gene,pw_g_list[[i]])
    }
    all_gene = unique(all_gene)

    useful_genes = intersect(all_gene,rownames(exp_mat))
    exp_mat = exp_mat[useful_genes,] # 选择pathway基因
    exp_mat = na.omit(exp_mat)

    gsva_matrix <- gsva(as.matrix(exp_mat), 
                        min.sz=min.sz, max.sz=max.sz,
                        pw_g_list,method='ssgsea', 
                        kcdf='Gaussian', 
                        abs.ranking=TRUE, parallel.sz=1)
    return(gsva_matrix)
}


# demo
if(FALSE){

}



