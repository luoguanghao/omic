#
#
get_depmap_subtype_sample <- function(info_df,mut_df=NULL,exp_df=NULL,subtype_term='lineage_subtype',subtype_name='AML'){
    sp_sample = (info_df[info_df[[subtype_term]]==subtype_name,])$DepMap_ID
    if(is.null(mut_df)&&is.null(exp_df)){
        return( na.omit(sp_sample) )
    }else if (!is.null(mut_df)) {
       return(mut_df%>%filter(DepMap_ID%in%sp_sample))
    }else if (!is.null(exp_df)) {
       return(exp_df[,c(colnames(exp_df)[1], 
                        intersect(colnames(exp_df)[2:ncol(exp_df)],
                                sp_sample))])
    }
}



if(FALSE){
    library(tidyverse)
    #########################################
    # import cleaned data
    #########################################
    mut_df = read_tsv('~/my_project/data/depmap/cleaned/CCLE_mutations_cleaned.tsv')
    exp_df = read_tsv('~/my_project/data/depmap/cleaned/CCLE_expression_cleaned.tsv')
    crispr_df = read_tsv('~/my_project/data/depmap/cleaned/CRISPR_gene_effect_cleaned.tsv')
    info_df = read_tsv('~/my_project/data/depmap/cleaned/sample_info.tsv')

    exp_mat = exp_df%>%column_to_rownames('Hugo_Symbol')

    #########################################
    #########################################

    ##
    # 为depmap多组学数据的整合设计一个模式
    # 所有'-'变成点'.'

    # get depmap data #########################################

    ## sample info
    info_df = read.csv('/y/Archive/Bondi/data/depmap/info/sample_info.csv',sep=',')



    ## exp mat

    exp_df = read_csv('~/my_project/data/depmap/expression/CCLE_expression.csv')
    colnames(exp_df)[1] = 'DepMap_ID'
    colnames(exp_df) = sapply(strsplit(colnames(exp_df),' \\('),'[[',1)
    exp_df = data.frame(t(exp_df%>%column_to_rownames('DepMap_ID')))%>%rownames_to_column('Hugo_Symbol')
    rownames(exp_df) = NULL
    exp_mat = exp_df%>%column_to_rownames('Hugo_Symbol')

    ## mut



    ## drug sen # 参照vizome格式 # 待完善
    drug_sen_df = read.csv('/y/Archive/Bondi/data/depmap/drug_sensitivity/sanger-dose-response.csv')
    gdsc_cell_ls = drug_sen_df[['Cell.line.name']]

    ## metabolics



    ## crispr
    # 读入depmap的crispr文件，给出一个矩阵，每一列是一个基因，每一行是样本
    #
    crispr_path = '~/my_project/data/depmap/dependence/CRISPR_gene_effect.csv'
    crispr_df = read_csv(crispr_path)
    colnames(crispr_df) = sapply(strsplit(colnames(crispr_df),' \\('),'[[',1)
    crispr_df$DepMap_ID = str_replace(crispr_df$DepMap_ID,'-','.')

    # meta info


    ## get sp cell line
    sp_cell_line = info_df[ grep('AML',info_df[['Subtype']]), ][['DepMap_ID']]
    cb_exp_mat = merge(info_df['DepMap_ID'],exp_mat,by="DepMap_ID",all=F) %>% column_to_rownames('DepMap_ID') %>% t()
    cb_exp_mat = cb_exp_mat[,intersect(sp_cell_line,colnames(cb_exp_mat))] # cell line不一定有转录组信息，选出有转录组信息的

}






























