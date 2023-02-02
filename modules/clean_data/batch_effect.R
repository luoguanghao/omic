

# data_sets: list for exp_mat: gene*sample,frist col is Hugo_Symbol,unique Hugo_Symbol
#
merge_data_set <- function(data_sets){
    # library(DESeq2)
    cmb_exp_mat = data_sets[[1]]
    condition_ls = rep(names(data_sets)[1], dim(data_sets[[1]])[2]-1)
    for(i in 2:length(data_sets)){
        cmb_exp_mat = 
            merge(cmb_exp_mat,data_sets[[i]],by.x='Hugo_Symbol',by.y='Hugo_Symbol')
        condition_ls = c( condition_ls, rep(names(data_sets)[i], dim(data_sets[[i]])[2]-1) )
    }
    cmb_exp_mat = cmb_exp_mat%>%column_to_rownames('Hugo_Symbol')
    # c(rep('depmap',dim(depmap_exp_count)[2]-1),rep('vizome',dim(vizome_exp_count)[2]-1))
    # return(list(cmb_exp_mat,condition_ls))
    coldata = data.frame(condition= condition_ls, 
                        row.names = colnames(cmb_exp_mat))

    # ===

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = cmb_exp_mat,
                    colData = coldata,
                    design = ~ condition)

    dds <- DESeq2::DESeq(dds)
    DESeq2::resultsNames(dds)
    res <- DESeq2::results(dds)

    vsd <- DESeq2::vst(dds, blind=FALSE)

    pcaData <- DESeq2::plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    p=ggplot2::ggplot(pcaData, aes(scale(PC1), scale(PC2), color=condition))+
        ggplot2::geom_point(size=3)+
        #xlim(-12, 12) +
        #ylim(-10, 10) +
        ggplot2::xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ggplot2::ylab(paste0("PC2: ",percentVar[2],"% variance")) #+
        geom_text(aes(label=name),vjust=2)
    plot(p)
    #ggsave("myPCAWithBatchEffect.png")
    
    # ===

    vsd <- DESeq2::vst(dds, blind=FALSE)
    SummarizedExperiment::assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$condition)

    pcaData <- DESeq2::plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    p=ggplot2::ggplot(pcaData, aes(PC1, PC2, color=condition)) +
        ggplot2::geom_point(size=3) +
        #xlim(-12, 12) +
        #ylim(-10, 10) +
        ggplot2::xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ggplot2::ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        ggplot2::geom_text(aes(label=name),vjust=2)
    plot(p)

    res_df = as.data.frame(res)

    normalized_exp_mat = SummarizedExperiment::assay(vsd)
    normalized_exp_mat = as.data.frame(normalized_exp_mat)

    return(list(normalized_exp_mat=normalized_exp_mat,
                deseq2_res_df=res_df,
                condition_ls=condition_ls,
                cmb_exp_mat=cmb_exp_mat
                ))
}














