

if(FALSE){
    exp_df = read_tsv('~/my_project/data/depmap/cleaned/CCLE_expression_cleaned.tsv')
    crispr_df = read_tsv('~/my_project/data/depmap/cleaned/CRISPR_gene_effect_cleaned.tsv')
    info_df = read_tsv('~/my_project/data/depmap/cleaned/sample_info.tsv')
    exp_mat = exp_df%>%column_to_rownames('Hugo_Symbol')

    cancer_ls = c('breast_carcinoma','colorectal_adenocarcinoma','neuroblastoma','bladder_carcinoma','ALL')
}
# 这里，根据某基因的CRISPR ceres分数，对样本进行分组，然后组间比较表达差异；可以限定在某个肿瘤背景下进行分析
# ！！！统一 `group_df`
#
# crispr_df:要求第一列是 DepMap_ID
# sel_crispr_df 可以外部给
do_compare_group_by_depmap_ceres <- function(subtype='pancancer', exp_mat, info_df, 
            crispr_df, mut_df, crispr_gene, th_l_1=-Inf,th_l_2=-1.5,th_r_1=-0.5,th_r_2=Inf,
            out_dir=NULL, pj_name=NULL,sel_crispr_df=NULL){
    source('/home/lgh/my_project/omic/modules/bioinformatics_wrapper/EDA.R')
    library(maftools)
    
    if(is.null(sel_crispr_df)){
        # 这里改一改能够根据药敏分组比较
        if(subtype=='pancancer'){    
            sel_crispr_df = (
                    crispr_df%>%filter(DepMap_ID%in%(info_df)$DepMap_ID)
                )[,c('DepMap_ID',crispr_gene)]
        }else{
            sel_crispr_df = (
                    crispr_df%>%filter(DepMap_ID%in%(info_df%>%filter(lineage_subtype==subtype))$DepMap_ID)
                )[,c('DepMap_ID',crispr_gene)]
        }

        sel_crispr_df$group = '0'

        sel_crispr_df[(th_l_1 < sel_crispr_df[[crispr_gene]])&(sel_crispr_df[[crispr_gene]]< th_l_2),
                    ]$group = sprintf('[%s,%s]',th_l_1,th_l_2)
        sel_crispr_df[(th_r_1 < sel_crispr_df[[crispr_gene]])&(sel_crispr_df[[crispr_gene]]< th_r_2),
                    ]$group = sprintf('[%s,%s]',th_r_1,th_r_2)

        sel_crispr_df = sel_crispr_df[order(sel_crispr_df[[crispr_gene]]),]
    }else{
        print(dim(sel_crispr_df))
    }


    # 表达差异分析
    condition_table = (sel_crispr_df%>%filter(group!='0'&DepMap_ID%in%colnames(exp_mat)))$group
    mat_for_limma = exp_mat[,(sel_crispr_df%>%filter(group!='0'&DepMap_ID%in%colnames(exp_mat)))$DepMap_ID]
    limma_res = do_limma(mat_for_limma, condition_table)
    
    # return(limma_res)

    if(is.null(out_dir)==FALSE){
        write_tsv(sel_crispr_df, 
                file=sprintf('%s/%s_crispr_score.tsv',out_dir,paste(c(subtype,pj_name),collapse='_'))
                )
        write_tsv(limma_res$resdata , 
                file=sprintf('%s/%s_limma_res.tsv',out_dir,paste(c(subtype,pj_name),collapse='_'))
                )
    }

    # 突变差异分析

    group_df = sel_crispr_df[,c(1,3)]%>%filter(group!='0')
    colnames(group_df) = c('sample','group')
    mut.df = mut_df%>%filter(DepMap_ID%in%group_df$sample)    
    colnames(group_df)[1] = 'Tumor_Sample_Barcode'

    maf = read.maf(maf = mut.df, clinicalData = group_df[,c(1,2)],
                vc_nonSyn=names(tail(sort(table(mut.df$Variant_Classification )))),
                verbose = FALSE)    

    oncoplot(
        maf = maf,
        clinicalFeatures = 'group',
        sortByAnnotation = TRUE)

    fab.ce = clinicalEnrichment(maf = maf, clinicalFeature = 'group')

    if(is.null(out_dir)==FALSE){
        write_tsv(fab.ce$pairwise_comparision, 
                file=sprintf('%s/%s_mut_pairwise_comparision.tsv',out_dir,paste(c(subtype,pj_name),collapse='_'))
                )
    }


    return(list(
        group_df=group_df,
        sel_crispr_df=sel_crispr_df,
        limma=list(limma_res=limma_res,
            mat_for_limma=mat_for_limma,
            condition_table=condition_table),
        mutdiff=list(mut.df=mut.df,
            pairwise_comparision=fab.ce$pairwise_comparision,
            groupwise_comparision=fab.ce$groupwise_comparision,
            fab.ce=fab.ce)
    ))
}


if(FALSE){
    
}





