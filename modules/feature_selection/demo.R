


# 1
source('/home/lgh/my_project/omic/modules/feature_selection/by_var.R')
source('/home/lgh/my_project/omic/modules/feature_selection/cor_one_by_one.R')
drugsen_gene_res = drugsen_gene_cor(drug_sen_df,
    exp_df,drug='ABT-263 (Navitoclax)',
    genes=get_most_var_feature(
        exp_df%>%column_to_rownames('Hugo_Symbol'),6000))







