
## do two group gsea


## do pre rank gsea





## 制作gsea需要的文件
# normalized_exp_mat
# condition_table
# comp_name 这一次工作的名字
prepare_gsea_file <- function(normalized_exp_mat,condition_table,result_dir,comp_name){
    # print(sprintf('## %s/%s_exp_mat.gct',result_dir,comp_name))
    gct_df = data.frame(NAME=rownames(normalized_exp_mat),
                        DESCRIPTION=rownames(normalized_exp_mat),
                        row.names=rownames(normalized_exp_mat))

    gct_df = cbind(gct_df,normalized_exp_mat)

    rownames(gct_df) = NULL

    gene_n = dim(gct_df)[1]
    sample_n = dim(gct_df)[2] - 2

    gct_header = list()
    gct_header[[colnames(gct_df)[1]]] = c('#',gene_n,colnames(gct_df)[1])
    gct_header[[colnames(gct_df)[2]]] = c('1.2',sample_n,colnames(gct_df)[2])

    for(ih in 3:dim(gct_df)[2]){
        gct_header[[colnames(gct_df)[ih]]] = c('','',colnames(gct_df)[ih])
    }
    gct_header = data.frame(gct_header)
    # return(list(gct_header,gct_df))
    
    gct_df = rbind(gct_header,gct_df)
    write_tsv(gct_df, file=sprintf('%s/%s_exp_mat.gct',result_dir,comp_name), col_names=FALSE)


    outfcon <- file(sprintf("%s/%s.cls",result_dir,comp_name), open="wt")
    outlines = paste(length(condition_table),length(unique(condition_table)),1,sep='\t')
    writeLines(outlines, con=outfcon)
    outlines = paste(c('#',unique(as.vector(condition_table))), collapse='\t')
    writeLines(outlines, con=outfcon) 
    outlines = paste(condition_table, collapse='\t')
    writeLines(outlines, con=outfcon) 
    close(outfcon) 
}


# pre rank GSEA：
# 提供rnk文件：
#
# # I will be ignore
# /gene name/ /foldChange/
# 
## need not to be sort!!!




# gct='mut/flt3_exp_mat.gct'
# cls='mut/flt3.cls'
# g1='wt'
# g2='mut'
# gmt='/home/lgh/my_project/data/pathway/h.all.v7.4.symbols.gmt'
# rpt_label='mut_flt_hallmark'
# metric='Signal2Noise'
# out='./mut/'
gsea_cmd <- function(gct, cls, g1, g2, gmt, rpt_label, metric, permute, out, scoring_scheme='weighted', num_show=100, win_linux='linux'){

    if(win_linux=='linux'){
        cmd = sprintf("/home/lgh/biosoft/GSEA_4.1.0/gsea-cli.sh GSEA \\
        -res %s \\
        -cls %s#%s_versus_%s \\
        -gmx %s \\
        -collapse No_Collapse \\
        -mode Max_probe \\
        -norm meandiv \\
        -nperm 1000 \\
        -permute %s \\
        -rnd_seed timestamp \\
        -rnd_type no_balance \\
        -scoring_scheme %s \\
        -rpt_label %s \\
        -metric %s \\
        -sort real \\
        -order descending \\
        -include_only_symbols true \\
        -make_sets true \\
        -median false \\
        -num 100 \\
        -plot_top_x %s \\
        -save_rnd_lists false \\
        -set_max 500 \\
        -set_min 15 \\
        -zip_report false \\
        -out %s",
        gct,cls,g2,g1,gmt,permute,scoring_scheme,rpt_label,metric,num_show,out)

        return(cmd)
    }else{

        cmd = sprintf("/home/lgh/biosoft/GSEA_4.1.0/gsea-cli.sh GSEA \\
        -res %s \\
        -cls %s#%s_versus_%s \\
        -gmx %s \\
        -collapse No_Collapse \\
        -mode Max_probe \\
        -norm meandiv \\
        -nperm 1000 \\
        -permute %s \\
        -rnd_seed timestamp \\
        -rnd_type no_balance \\
        -scoring_scheme %s \\
        -rpt_label %s \\
        -metric %s \\
        -sort real \\
        -order descending \\
        -include_only_symbols true \\
        -make_sets true \\
        -median false \\
        -num 100 \\
        -plot_top_x 50 \\
        -save_rnd_lists false \\
        -set_max 500 \\
        -set_min 15 \\
        -zip_report false \\
        -out %s",
        gct,cls,g2,g1,gmt,permute,scoring_scheme,rpt_label,metric,out)

        return(cmd)     

    }

}

gsea_cmd_prerank <- function(rnk, gmt, rpt_label, permute, out, scoring_scheme='classic', win_linux='linux'){
    cmd = sprintf("/home/lgh/biosoft/GSEA_4.1.0/gsea-cli.sh GSEAPreranked \\
        -gmx %s \\
        -collapse No_Collapse \\
        -mode Abs_max_of_probes \\
        -norm meandiv \\
        -nperm %s \\
        -rnd_seed timestamp \\
        -rnk %s \\
        -scoring_scheme %s \\
        -rpt_label %s \\
        -create_svgs false \\
        -include_only_symbols true \\
        -make_sets true \\
        -plot_top_x 20 \\
        -set_max 500 \\
        -set_min 15 \\
        -zip_report false \\
        -out %s",
    gmt,permute,rnk,scoring_scheme,rpt_label,out)

    return(cmd)
}


if(FALSE){
    ## demo run
    # from deseq2 result
    # 不能在mat里面有NA
    source('/home/lgh/my_project/omic/modules/bioinformatics_wrapper/prepare_gsea.R')
    source('/home/lgh/my_project/omic/project/hl/大组学/组学相关分析/cross_omics_analysis.R')

    group_df = data.frame(sample=colnames(depmap_exp_df)[-1], group='no_bc')
    group_df[group_df$sample%in%bc_cell_line_df$DepMap_ID,]$group = 'bc'

    normalized_exp_mat = bc_vs_other$resdata[,-c(1:7)]
    rownames(normalized_exp_mat) = bc_vs_other$resdata$row.names

    # gsea_res = list()
    #for(i in length(gmtfile)){
    #    term = names(gmtfile)[[i]]
    cmd = post_subtyping.prepare_gsea(normalized_exp_mat=normalized_exp_mat, group_df=group_df, 
                                g1='no_bc', g2='bc', 
                                out='res_0509/cell_line', 
                                gmtfile = gmtfile, 
                                result_dir='res_0509/cell_line', comp_name=sprintf('bc_vs_nobc_%s',term))
        


    gsea_res = list()
    for(i in 1:length(cmd)){
        gsea_res[[term]] = system(cmd[[i]], intern = TRUE)
    }







}