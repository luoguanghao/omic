


### for rnaseq from deseq2 to gsea
if(FALSE){

compare_list = list(
    c('NC_siC','NC_siD'),
    c('NC_siC','TGFb_siC'),
    c('TGFb_siC','TGFb_siD')
)
count_df = read_tsv('./feature_count_matrix.tsv')
coldata = data.frame(condition=c(rep('NC_siC',3),rep('NC_siD',3),rep('TGFb_siC',3),rep('TGFb_siD',3)),row.names=colnames(count_df)[-1])
group_df = data.frame(sample=colnames(count_df)[-1],group=c(rep('NC_siC',3),rep('NC_siD',3),rep('TGFb_siC',3),rep('TGFb_siD',3)))
gsea_result_dir = '/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae51/server/raw_data/lihao/analysis/new_gsea_result'

}

if(FALSE){

do_deseq2_new_multi_res = do_deseq2_new_multi(count_df%>%column_to_rownames('Hugo_Symbol'), group_df )


for(comp in compare_list){

    deg_df = data.frame(results(do_deseq2_new_multi_res$dds,contrast=c("condition",comp[1],comp[2])))
    
    normalized_exp_mat = log2(1+counts(do_deseq2_new_multi_res$dds,normalize=TRUE))
    normalized_exp_mat = normalized_exp_mat[rownames(deg_df%>%filter(baseMean>=10)),]

    # group_df = coldata%>%rownames_to_column('sample')
    colnames(group_df)[2] = 'group'

    tmp_group_df = group_df%>%filter(group%in%comp)

    post_subtyping.prepare_gsea_res = post_subtyping.prepare_gsea(normalized_exp_mat[,(tmp_group_df)$sample], tmp_group_df, comp[2], comp[1], 
                                gsea_result_dir, gmtfile, gsea_result_dir, paste(compare_list[[1]], collapse='_vs_'), 100)

    for(cmd in post_subtyping.prepare_gsea_res){

        system(cmd)

    }

}



}

