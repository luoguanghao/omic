library(DESeq2)



##################
# 检查数据是否都合法
##################





# 包括各种转换，tpm转rpkm等等
# https://zhuanlan.zhihu.com/p/150300801
#
#
#############
# count2tpm
#############
countToTpm <- function(counts, effLen)
{
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
}



###########
# rpkm转tpm
###########
# raw_exp_mat: gene*sample rpkm without log
## out:
# exp_mat_tpm gene*sample, tpm without log, sum of every columns is 1e6
# ===
if(FALSE){
    fpkmToTpm <- function(fpkm)
    {
        exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
    }
    exp_mat_tpm = apply(raw_exp_mat,2,fpkmToTpm)
}










###############
### do vrt normalization with use of deseq2 
###############
# do vrt normalization with use of deseq2 #
# 
# deseq2_exp_mat: matrix ,rowname is gene, colnames is sample, must be integer
# condition_table: a vector define the group of every sample
# 
# possible process for deseq2
##
# condition_table = factor(c(rep('sen',5),rep('unsen',5)))
# deseq2_exp_mat = as.matrix(exp_mat%>%column_to_rownames('gene_id'))
# deseq2_exp_mat = apply(deseq2_exp_mat,2,as.integer)
# rownames(deseq2_exp_mat) = exp_mat$gene_id  // rownames(deseq2_exp_mat) = rownames(exp_mat)
deseq2_normalization <- function(deseq2_exp_mat,condition_table=NULL, to_int=FALSE){
    samples_num = length(colnames(deseq2_exp_mat))

    if(to_int){
        rn = rownames(deseq2_exp_mat)
        deseq2_exp_mat = sapply(deseq2_exp_mat, round)
        rownames(deseq2_exp_mat) = rn
    }

    if(is.null(condition_table)){
        condition_table = factor(c(rep('sen',round(samples_num/2)), rep('unsen',samples_num-round(samples_num/2))))
    }
    # =========
    dds <- DESeqDataSetFromMatrix(deseq2_exp_mat, 
                                DataFrame(condition_table), 
                                design= ~ condition_table)
    dds <- dds[rowSums(counts(dds)) > 1,]
    dds2 <- DESeq(dds)
    resultsNames(dds2)
    # acquire the results using function results(), and assign to res
    res <- results(dds2)
    # view the summary of results
    summary(res)
    # resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
    # output
    #head(resdata)

    deseq2_normalized_exp_mat = as.data.frame(counts(dds2,normalize=TRUE))
    vsd <- vst(dds2, blind=FALSE)

    return(list(norm_exp_mat=deseq2_normalized_exp_mat, 
                vst_exp_mat=assay(vsd),res=res, dds=dds2))
}













