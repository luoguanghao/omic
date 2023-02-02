


#########################################
# get TCGA data #########################
#########################################

## exp mat #########################################
if(FALSE){
    exp_mat = read.csv('/y/Archive/Bondi/data/TCGA/LAML/TCGA-LAML.htseq_fpkm.tsv.gz',sep='\t')
    probeMap_df = read.csv('/y/Archive/Bondi/data/TCGA/LAML/gencode.v22.annotation.gene.probeMap',sep='\t')
}
###

load_tcga_exp_mat <- function(exp_mat, probeMap_df){
    exp_mat = read_tsv(exp_path)
    probeMap_df = read.csv(probeMap_path,sep='\t')
    ###
    probeMap_df <- probeMap_df %>%  dplyr::select(1,2)
    colnames(probeMap_df) <- c("Ensembl_ID","gene_name")

    exp_mat[,c(2:dim(exp_mat)[2])] = 2**exp_mat[-1]-1 # for count data, reverse

    exp_mat <- inner_join(probeMap_df, exp_mat, by = "Ensembl_ID") %>% dplyr::select(-1)
    exp_mat <- aggregate(.~ gene_name, exp_mat, mean)
    rn = rownames(exp_mat)
    exp_mat <- exp_mat %>% column_to_rownames("gene_name")
    exp_mat = apply(exp_mat,2,as.integer)
    row.names(exp_mat) = rn
    return(exp_mat)
}

