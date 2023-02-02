

get_most_var_feature <- function(mat,top_n=5000,genes_vec='genes'){
    feature_ls = apply(mat,1,var)
    feature_ls = feature_ls[order(feature_ls,decreasing=TRUE)]
    if(genes_vec!='genes') return(feature_ls)
    return(names(feature_ls[1:top_n]))
}





