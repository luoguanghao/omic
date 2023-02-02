get_time_suffix <- function(){
    return(str_replace_all(str_replace_all(str_replace_all(Sys.time(),'-',''),' ','_'),':',''))
}

if(FALSE){ # resource
    # 组学数据

    # signature gene set
}

################################
# pipeline example:vizome
################################

## 代谢分型
if(FALSE){
    library(tidyverse)
    library(DESeq2)
    library(GSVA)
    library(pheatmap)
    library(ConsensusClusterPlus)
    library(ggfortify)

    # source('/home/lgh/my_project/omic/modules/prepare_data/depmap_data.R')
    # source('/home/lgh/my_project/omic/project/protostatsis/func_subtyping.R')
    source('/home/lgh/my_project/omic/modules/grouping/func_subtyping.R')
    source('/home/lgh/my_project/omic/modules/bioinformatics_wrapper/gsva.R')
    source('/home/lgh/my_project/omic/modules/grouping/grouping_by_clustering_pipeline.R')
    source('/home/lgh/my_project/omic/modules/grouping/func_post_subtyping.R')
    # =====================

    project_name = 'vizome_metabolism_subtyping'
    result_dir = sprintf('./result_test_pipe_line_vizome_%s/',get_time_suffix())
    dir.create(result_dir)
    # subtyping_res = list()


    ## info data
    #info_df_path = '~/my_project/data/depmap/cleaned/sample_info.tsv'
    #info_df = read_tsv(info_df_path)
    ## exp data
    #exp_df_path = '~/my_project/data/depmap/cleaned/CCLE_expression_cleaned.tsv'
    #exp_df = read_tsv(exp_df_path)
    #normalized_exp_mat = get_depmap_subtype_sample(info_df=info_df,exp_df=exp_df,subtype_term='lineage_subtype',subtype_name='AML')%>%
    #                        column_to_rownames('Hugo_Symbol')

    ## exp data
    exp_df = read_tsv('~/my_project/data/VIZOME/cleaned/VIZOME_log2fpkm_expression_cleaned.tsv')
    normalized_exp_mat = exp_df%>%column_to_rownames('Hugo_Symbol')



    # pathway info 
    pathway_ls_path = './data/MetabolicPathway_ls.txt'
    pw_data <- read.csv(pathway_ls_path, encoding="UTF-8",skip=0,sep='\t',stringsAsFactors = FALSE)
    ### 

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

    ####
    # go
    ####

    gsva_matrix = do_gsva(exp_mat=normalized_exp_mat,pw_data=pw_data,pw_g_list=pw_g_list)

    subtyping_res = do_ConsensusCluster(exp_mat_cc=gsva_matrix,
        maxK_sample=8,maxK_gene=8,reps=10000,
        result_dir=NULL,project_name=NULL)

    pathway_anno_df = read_tsv('./data/pathway_anno_df.txt')
    pathway_anno_df = pathway_anno_df%>%column_to_rownames('pathway')
    colnames(pathway_anno_df)[1] = 'Cluster_row'

    subtyping_res = do_ConsensusCluster_result(subtyping_res=subtyping_res,
        n_cluster_sample=4, n_cluster_gene=4,
        rownames_define=NULL,colnames_define=NULL,
        annotation_row_define=pathway_anno_df, annotation_col_define=NULL,
        result_dir=result_dir, project_name=project_name,hmap_label='pathway_anno')

    # ========================
    # post subtyping
    # ========================
    suv_df = read_tsv('~/my_project/data/VIZOME/cleaned/VIZOME_survival_cleaned.tsv')

    gmtfile = list(
        go=paste(PATH_DATA,'/pathway/c5.go.bp.v7.2.symbols.gmt',sep=''),
        kegg=paste(PATH_DATA,'/pathway/c2.cp.kegg.v7.2.symbols.gmt',sep=''),
        reactome=paste(PATH_DATA,'/pathway/reactome.gmt',sep='')
    )

    maf.df = read_tsv('/home/lgh/my_project/data/VIZOME/cleaned/VIZOME_mutation.maf')

    drug_family_df = read_tsv('~/my_project/data/VIZOME/cleaned/VIZOME_drug_families.tsv')
    drug_anno_df = read_tsv('~/my_project/data/VIZOME/cleaned/VIZOME_drug_annotation.tsv')
    drug_sen_df = read_tsv('~/my_project/data/VIZOME/cleaned/VIZOME_drug_sen_cleaned.tsv')

    suv_data_list=list(suv_df=suv_df)
    exp_data_list=list(exp_mat=normalized_exp_mat,gmtfile=gmtfile)
    mut_data_list=list(maf.df=maf.df)
    drug_data_list=list(drug_sen_df=drug_sen_df,drug_anno_df=drug_anno_df,drug_family_df=drug_family_df,ic50auc='auc')
    other_subtyping_data_list=list()

    # 流程分析：生存，突变，表达，药敏，不同分型的比较
    #  group_df
    res_post_subtyping = post_subtyping(subtyping_res=subtyping_res,group_df=NULL,
            suv_data_list=suv_data_list,exp_data_list=exp_data_list,
            mut_data_list=mut_data_list,drug_data_list=drug_data_list,
            other_subtyping_data_list=NULL)

    #############
    # 个性化分析 #
    #############

    ## 突变的个性化分析 ##
    
    # 突变综合分析
    compare_maf_res = compare_maf(mut_data_list$maf.df, subtyping_res$sample_group_df)

    # 突变通路分析

    # 突变pair分析
    pair_compare_maf_res = pair_compare_maf(maf.df=mut_data_list$maf.df, group_df=subtyping_res$sample_group_df,
            m1Name='c1', m2Name='c3', m1show='c1', m2show='c3',
            my_genes=NULL, sig_th=0.05, plot=TRUE)

    ## 药敏的个性化分析 ##

    # 两个group进行药敏比较比较
    compare_drug_sen_two_drug_res = compare_drug_sen_two_drug(comp_res=res_post_subtyping$drug$comp_res,
                            drug_family_df=drug_data_list$drug_family_df,
                            g1='c1',g2='c3',
                            g1_name='g1 carbohydrate sen',g2_name='g3 lipoid sen')
    compare_drug_sen_two_drug_res
    # 药物靶点的富集
    ## GSEA类型的富集
    ## ORA类型的富集



}

################################
# pipeline example depmap
################################


if(FALSE){
    library(tidyverse)
    library(DESeq2)
    library(GSVA)
    library(pheatmap)
    library(ConsensusClusterPlus)
    library(ggfortify)

    # source('/home/lgh/my_project/omic/modules/prepare_data/depmap_data.R')
    source('/home/lgh/my_project/omic/project/protostatsis/func_subtyping.R')
    source('/home/lgh/my_project/omic/modules/grouping/func_subtyping.R')
    source('/home/lgh/my_project/omic/modules/bioinformatics_wrapper/gsva.R')
    source('/home/lgh/my_project/omic/modules/grouping/grouping_by_clustering_pipeline.R')

    source('/home/lgh/my_project/omic/modules/prepare_data/depmap_data.R')

    project_name = 'depmap_metabolism_subtyping'
    result_dir = sprintf('./result_test_pipe_line_depmap_%s/',get_time_suffix())
    dir.create(result_dir)
    # subtyping_res = list()

    info_df = read_tsv('/home/lgh/my_project/data/depmap/cleaned/sample_info.tsv')

    ## exp data
    exp_df = read_tsv('/home/lgh/my_project/data/depmap/cleaned/CCLE_expression_cleaned.tsv')
    exp_df = get_depmap_subtype_sample(info_df,mut_df=NULL,exp_df=exp_df,subtype_term='lineage_subtype',subtype_name='AML')
    normalized_exp_mat = exp_df%>%column_to_rownames('Hugo_Symbol')



    # pathway info 
    pathway_ls_path = './data/MetabolicPathway_ls.txt'
    pw_data <- read.csv(pathway_ls_path, encoding="UTF-8",skip=0,sep='\t',stringsAsFactors = FALSE)
    ### 

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

    ####
    # go
    ####

    gsva_matrix = do_gsva(exp_mat=normalized_exp_mat,pw_data=pw_data,pw_g_list=pw_g_list)

    subtyping_res = do_ConsensusCluster(exp_mat_cc=gsva_matrix,
        maxK_sample=8,maxK_gene=8,reps=10000,
        result_dir=NULL,project_name=NULL)

    pathway_anno_df = read_tsv('./data/pathway_anno_df.txt')
    pathway_anno_df = pathway_anno_df%>%column_to_rownames('pathway')
    colnames(pathway_anno_df)[1] = 'Cluster_row'

    define_col_vec = (info_df[,c('DepMap_ID','stripped_cell_line_name')]%>%column_to_rownames('DepMap_ID'))[colnames(subtyping_res$exp_mat),]
    names(define_col_vec) = colnames(subtyping_res$exp_mat)

    subtyping_res = do_ConsensusCluster_result(subtyping_res=subtyping_res,
        n_cluster_sample=4, n_cluster_gene=4,
        rownames_define=NULL,colnames_define=define_row_col_vec,
        annotation_row_define=pathway_anno_df, annotation_col_define=NULL,
        result_dir=result_dir, project_name=project_name,hmap_label='pathway_anno')

    # ========================
    # post subtyping
    # ========================


}



################################
# pipeline example depmap
################################


if(FALSE){
    library(tidyverse)
    library(DESeq2)
    library(GSVA)
    library(pheatmap)
    library(ConsensusClusterPlus)
    library(ggfortify)

    # source('/home/lgh/my_project/omic/modules/prepare_data/depmap_data.R')
    source('/home/lgh/my_project/omic/project/protostatsis/func_subtyping.R')
    source('/home/lgh/my_project/omic/modules/grouping/func_subtyping.R')
    source('/home/lgh/my_project/omic/modules/bioinformatics_wrapper/gsva.R')
    source('/home/lgh/my_project/omic/modules/grouping/grouping_by_clustering_pipeline.R')

    source('/home/lgh/my_project/omic/modules/prepare_data/depmap_data.R')

    project_name = 'depmap_metabolism_subtyping'
    result_dir = sprintf('./result_test_pipe_line_depmap_%s/',get_time_suffix())
    dir.create(result_dir)
    # subtyping_res = list()

    info_df = read_tsv('/home/lgh/my_project/data/depmap/cleaned/sample_info.tsv')

    ## exp data
    exp_df = read_tsv('/home/lgh/my_project/data/depmap/cleaned/CCLE_expression_cleaned.tsv')
    exp_df = get_depmap_subtype_sample(info_df,mut_df=NULL,exp_df=exp_df,subtype_term='lineage_subtype',subtype_name='AML')
    normalized_exp_mat = exp_df%>%column_to_rownames('Hugo_Symbol')



    # pathway info 
    pathway_ls_path = './data/MetabolicPathway_ls.txt'
    pw_data <- read.csv(pathway_ls_path, encoding="UTF-8",skip=0,sep='\t',stringsAsFactors = FALSE)
    ### 

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

    ####
    # go
    ####

    gsva_matrix = do_gsva(exp_mat=normalized_exp_mat,pw_data=pw_data,pw_g_list=pw_g_list)

    subtyping_res = do_ConsensusCluster(exp_mat_cc=gsva_matrix,
        maxK_sample=8,maxK_gene=8,reps=10000,
        result_dir=NULL,project_name=NULL)

    pathway_anno_df = read_tsv('./data/pathway_anno_df.txt')
    pathway_anno_df = pathway_anno_df%>%column_to_rownames('pathway')
    colnames(pathway_anno_df)[1] = 'Cluster_row'

    define_col_vec = (info_df[,c('DepMap_ID','stripped_cell_line_name')]%>%column_to_rownames('DepMap_ID'))[colnames(subtyping_res$exp_mat),]
    names(define_col_vec) = colnames(subtyping_res$exp_mat)

    subtyping_res = do_ConsensusCluster_result(subtyping_res=subtyping_res,
        n_cluster_sample=4, n_cluster_gene=4,
        rownames_define=NULL,colnames_define=define_row_col_vec,
        annotation_row_define=pathway_anno_df, annotation_col_define=NULL,
        result_dir=result_dir, project_name=project_name,hmap_label='pathway_anno')

    # ========================
    # post subtyping
    # ========================

    # 

    
}

################################
# pipeline example ch
################################


if(FALSE){
    library(tidyverse)
    library(DESeq2)
    library(GSVA)
    library(pheatmap)
    library(ConsensusClusterPlus)
    library(ggfortify)

    # source('/home/lgh/my_project/omic/modules/prepare_data/depmap_data.R')
    source('/home/lgh/my_project/omic/project/protostatsis/func_subtyping.R')
    source('/home/lgh/my_project/omic/modules/grouping/func_subtyping.R')
    source('/home/lgh/my_project/omic/modules/bioinformatics_wrapper/gsva.R')
    source('/home/lgh/my_project/omic/modules/grouping/grouping_by_clustering_pipeline.R')

    source('/home/lgh/my_project/omic/modules/prepare_data/depmap_data.R')

    project_name = 'depmap_metabolism_subtyping'
    result_dir = sprintf('./result_test_pipe_line_depmap_%s/',get_time_suffix())
    dir.create(result_dir)
    # subtyping_res = list()

    info_df = read_tsv('/home/lgh/my_project/data/depmap/cleaned/sample_info.tsv')

    ## exp data
    exp_df = read_tsv('/home/lgh/my_project/data/depmap/cleaned/CCLE_expression_cleaned.tsv')
    exp_df = get_depmap_subtype_sample(info_df,mut_df=NULL,exp_df=exp_df,subtype_term='lineage_subtype',subtype_name='AML')
    normalized_exp_mat = exp_df%>%column_to_rownames('Hugo_Symbol')



    # pathway info 
    pathway_ls_path = './data/MetabolicPathway_ls.txt'
    pw_data <- read.csv(pathway_ls_path, encoding="UTF-8",skip=0,sep='\t',stringsAsFactors = FALSE)
    ### 

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

    ####
    # go
    ####

    gsva_matrix = do_gsva(exp_mat=normalized_exp_mat,pw_data=pw_data,pw_g_list=pw_g_list)

    subtyping_res = do_ConsensusCluster(exp_mat_cc=gsva_matrix,
        maxK_sample=8,maxK_gene=8,reps=10000,
        result_dir=NULL,project_name=NULL)

    pathway_anno_df = read_tsv('./data/pathway_anno_df.txt')
    pathway_anno_df = pathway_anno_df%>%column_to_rownames('pathway')
    colnames(pathway_anno_df)[1] = 'Cluster_row'

    define_col_vec = (info_df[,c('DepMap_ID','stripped_cell_line_name')]%>%column_to_rownames('DepMap_ID'))[colnames(subtyping_res$exp_mat),]
    names(define_col_vec) = colnames(subtyping_res$exp_mat)

    subtyping_res = do_ConsensusCluster_result(subtyping_res=subtyping_res,
        n_cluster_sample=4, n_cluster_gene=4,
        rownames_define=NULL,colnames_define=define_row_col_vec,
        annotation_row_define=pathway_anno_df, annotation_col_define=NULL,
        result_dir=result_dir, project_name=project_name,hmap_label='pathway_anno')

    # ========================
    # post subtyping
    # ========================


}







## 普通基因分型
if(FALSE){

}


## xbp1分型
if(FALSE){

}







