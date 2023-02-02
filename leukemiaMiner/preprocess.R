

######################
# select sample by NA
######################
# pro_mat_raw,pro_mat_impute: sample*gene or sample*gene
# return_type: genes or matrix
# s_g: gene*sample or sample*gene
# precent of NA be saved
select_gene_according_to_NA <- function(pro_mat_raw,pro_mat_impute,precent=0.5, s_g='gs',return_type='genes'){

    sel_pro_by_na = rownames(pro_mat_raw)[apply(pro_mat_raw,1,function(x) sum(is.na(x)))<ncol(pro_mat_raw)*precent]

    if(return_type=='genes'){
        return(sel_pro_by_na)
    }
}



######################
# deal NA
######################
# MOVICS有这个功能，但是不够用
# min，knn，补0法进行补miss










