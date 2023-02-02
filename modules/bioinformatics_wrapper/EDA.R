# limma
# DESeq2
# egeR
# http://localhost:9424/notebooks/myhome/my_project/jf/rnaseq_9_15/jf_rnaseq.ipynb 有

# limma #
# 其规则是condition_table 中后比前获得fc值
# condition_table是一个vec,对应mat的样本,标识group
# mat_for_limma是mat: gene*sample
#
# result:
#
do_limma <- function(mat_for_limma, condition_table) {
    cat(sprintf('foldchange %s : %s',unique(as.vector(condition_table))[2],
                unique(as.vector(condition_table))[1]))

    library(limma)
    design <- model.matrix(~condition_table)
    colnames(design) <- levels(condition_table)
    rownames(design) <- colnames(mat_for_limma)
    fit <- lmFit(mat_for_limma, design)
    fit <- eBayes(fit, trend=TRUE)
    limma_res <- topTable(fit, coef=2,n=Inf)

    limma_res = cbind(limma_res,mat_for_limma[rownames(limma_res),])
    limma_res=limma_res%>%rownames_to_column('row.names')
    colnames(limma_res)[1:7] = c('row.names','log2FoldChange','AveExpr','stat','p','padj','B')
    
    return(list(resdata=limma_res))

}

# mat_for_limma是mat: gene*sample
# group_df: $sample ; $group 两列
# g2 vs g1 后比前
# FC是两组的mean值相减
# limma需要log的矩阵
do_limma_new <- function(mat_for_limma, group_df, g1, g2) {
    cat(sprintf('foldchange %s : %s\n',g2,g1))
    print(table(group_df$group))
    group_df = data.frame(group_df)
    group_df = group_df%>%filter(sample%in%colnames(mat_for_limma))

    group_df$group<-as.character(group_df$group)
    
    group_df = rbind(group_df%>%filter(group%in%c(g1)),group_df%>%filter(group%in%c(g2)))
    
    mat_for_limma = mat_for_limma[,group_df$sample]
    
    group_df[group_df$group==g2,]$group = 2
    group_df[group_df$group==g1,]$group = 1
    condition_table = group_df$group
    
    library(limma)
    design <- model.matrix(~condition_table)
    colnames(design) <- levels(condition_table)
    rownames(design) <- colnames(mat_for_limma)
    fit <- lmFit(mat_for_limma, design)
    fit <- eBayes(fit, trend=TRUE)
    limma_res <- topTable(fit, coef=2,n=Inf)

    limma_res = cbind(limma_res,mat_for_limma[rownames(limma_res),])
    limma_res=limma_res%>%rownames_to_column('row.names')
    colnames(limma_res)[1:7] = c('row.names','log2FoldChange','AveExpr','stat','p','padj','B')
    
    return(list(resdata=limma_res))

}


##################################################################
##################################################################

# DESeq2 #
# 其规则是condition_table 中后比前获得fc值
# group_df : sample;group 
# deseq2_exp_mat是mat: gene*sample
#
#- 输出矩阵，样本的顺序就是g1,g2的顺序
#- 是数字大的比数字小的 condition table中，所以这就是g2 vs g1
#
# result:
#

do_deseq2_new_multi <- function(deseq2_exp_mat, group_df, compare=NULL, to_int=FALSE) {
    coldata = group_df
    colnames(coldata)[2] = 'condition'
    coldata = coldata%>%column_to_rownames('sample')
    dds <- DESeqDataSetFromMatrix(countData = deseq2_exp_mat,
                                colData = coldata,
                                design = ~ condition)
    dds <- DESeq(dds)
    norm_df = data.frame(counts(dds,normalize=TRUE))%>%rownames_to_column('gene')
    vst_df = data.frame(SummarizedExperiment::assay(vst(dds, blind=FALSE)))%>%rownames_to_column('gene')

    # if(compare!=NULL)

    return(list(norm_df=norm_df, vst_df=vst_df, dds=dds))

}

do_deseq2_new <- function(deseq2_exp_mat, group_df, g1, g2, to_int=FALSE) {
    if(to_int==TRUE){
        rn = rownames(deseq2_exp_mat)
        deseq2_exp_mat = sapply(deseq2_exp_mat,round)
        rownames(deseq2_exp_mat) = rn
    }

    cat(sprintf('foldchange %s : %s',g2,g1))

    group_df = rbind(group_df%>%filter(group%in%c(g1)),group_df%>%filter(group%in%c(g2)))
    group_df = group_df%>%filter(sample%in%colnames(deseq2_exp_mat))
    
    deseq2_exp_mat = deseq2_exp_mat[,group_df$sample]

    group_df[group_df$group==g1,]$group = '1'   
    group_df[group_df$group==g2,]$group = '2'

    condition_table = group_df$group
    
    library(DESeq2)
    # deseq2_exp_mat

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
    resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=TRUE)

    vsd <- vst(dds2, blind=FALSE)
    resdata_vst <- merge(as.data.frame(res),SummarizedExperiment::assay(vsd),by="row.names",sort=TRUE)

    # output
    resdata = resdata[order(resdata$padj),]
    resdata_vst = resdata_vst[order(resdata_vst$padj),]

    colnames(resdata)[1] = 'row.names'
    colnames(resdata_vst)[1] = 'row.names'

    return(list(resdata=resdata, resdata_vst=resdata_vst, dds=dds2))
}


do_deseq2 <- function(deseq2_exp_mat, condition_table) {
    cat(sprintf('foldchange %s : %s',unique(as.vector(condition_table))[2],
                unique(as.vector(condition_table))[1]))
    
    library(DESeq2)
    # deseq2_exp_mat

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
    resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=TRUE)

    vsd <- vst(dds2, blind=FALSE)
    resdata_vst <- merge(as.data.frame(res),SummarizedExperiment::assay(vsd),by="row.names",sort=TRUE)

    # output
    resdata = resdata[order(resdata$padj),]
    resdata_vst = resdata_vst[order(resdata_vst$padj),]

    colnames(resdata)[1] = 'row.names'
    colnames(resdata_vst)[1] = 'row.names'

    return(list(resdata=resdata, resdata_vst=resdata_vst, dds=dds2))
}

## 朴素版本的多组比较 deseq2
# deseq2_exp_mat是mat: gene*sample
# anno_df:描述每个样本的归属:group,sample_name两列
do_deseq2_multi_compare <- function(exp_mat, anno_df){
    comp = as.data.frame(combn(levels(anno_df$group),2))
    res_list = list()
    for(c in colnames(comp)){
        sel_sample = (anno_df%>%filter(group%in%comp[[c]]))$sample_name
        condition_table = (anno_df%>%filter(group%in%comp[[c]]))$group
        
        res_list[[paste(comp[[c]],collapse='_')]] = res = do_deseq2(deseq2_exp_mat=exp_mat[,sel_sample], condition_table)
    }

    return(res_list)
}


## deseq2多组比较
#
# deseq2_exp_mat是mat: gene*sample
# anno_df:描述每个样本的归属
# 
## 读取结果: name= 处填入resultsNames(dds2)的compare name
# res <- results(dds2, name="condition_jeko1sgctr0_vs_jeko1sgctrldm")
# resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=TRUE)
# head(resdata)
#
do_multi_deseq2 <- function(exp_mat,anno_df,ref){

    coldata <- data.frame(
                condition = factor(anno_df$group)
                )
    # set `ctl` as reference level
    coldata$condition <- relevel(coldata$condition, ref = ref)


    dds <- DESeqDataSetFromMatrix(
    countData = exp_mat[,anno_df$sample_name],
    colData = coldata,
    design= ~ condition)

    dds2 <- DESeq(dds, betaPrior = FALSE)
    print(resultsNames(dds2))

    return(dds2)
}




















