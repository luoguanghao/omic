source('/home/lgh/my_project/omic/modules/grouping/func_subtyping.R')

# discovery_result

discovery.CC.cluster.first_step <- function(exp_mat_cc,scale=TRUE,
        maxK_sample=8,maxK_gene=8,
        clusterAlg='km',distance='euclidean',
        pItem=0.8,reps=1000, pdf_png='pdf',
        result_dir=NULL,project_name=NULL){


    subtyping_res = do_ConsensusCluster(exp_mat_cc=normalized_exp_mat_for_cluster,
        maxK_sample=maxK_sample,maxK_gene=maxK_gene,reps=reps,
        result_dir=result_dir,project_name=project_name)

    discovery_result = list(
        method='cc_cluster',
        result=list(subtyping_res=subtyping_res),
        group_df=NULL
    )
    return(discovery_result)
}



discovery.CC.cluster.second_step <- function(discovery_result,
        n_cluster_sample=3, n_cluster_gene=3,
        rownames_define=NULL,colnames_define=NULL,
        annotation_row_define=NULL,annotation_col_define=NULL,
        result_dir=NULL,project_name=NULL,hmap_label=NULL,
        fontsize_row=10,fontsize_col=NULL,wid_height_scale=1,
        show_rownames=FALSE,show_colnames=FALSE,scale=TRUE){

    subtyping_res2 = do_ConsensusCluster_result(subtyping_res=discovery_result$result$subtyping_res,
        n_cluster_sample=n_cluster_sample, n_cluster_gene=n_cluster_gene,
        rownames_define=rownames_define,colnames_define=colnames_define,
        annotation_col_define=annotation_col_define,wid_height_scale = wid_height_scale,
        result_dir=result_dir, project_name=project_name)

    discovery_result = list(
        method='cc_cluster',
        result=list(subtyping_res=discovery_result$result$subtyping_res, subtyping_res2=subtyping_res2),
        group_df=subtyping_res2$sample_group_df
    )
    discovery_result$group_df$group = as.character(discovery_result$group_df$group)

    return(discovery_result)
}



discovery.HCPC.cluster <- function(exp_mat_cc,ncp=5, fig=FALSE){

    do_HCPC_cluster_res = do_HCPC_cluster(exp_mat_cc=exp_mat_cc, ncp=ncp, fig=fig)

    group_df = data.frame(sample=rownames(do_HCPC_cluster_res$res.hcpc$data.clust),
                            group=do_HCPC_cluster_res$res.hcpc$data.clust$clust)
    group_df$group = as.character(group_df$group)

    discovery_result = list(
        method='HCPC.cluster',
        result=do_HCPC_cluster_res,
        group_df=group_df
    )

    return(discovery_result)

}






















