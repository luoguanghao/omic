# 加上对deseq2的支持
#
# 输入 gmt文件，表达矩阵
# 输出 差异表达分析，富集分析结果表格，相应可视化；用于GSEA软件的表达矩阵和相关文件；用于网络分析的相关文件
#

library(tidyverse)
library(limma)

# limma #
# 其规则是condition_table 中后比前获得fc值
# condition_table是一个vec,对应mat的样本,标识group
# mat_for_limma是mat: gene*sample
do_limma <- function(mat_for_limma, condition_table) {
    design <- model.matrix(~condition_table)
    colnames(design) <- levels(condition_table)
    rownames(design) <- colnames(mat_for_limma)
    fit <- lmFit(mat_for_limma, design)
    fit <- eBayes(fit, trend=TRUE)
    limma_res <- topTable(fit, coef=2,n=Inf)
    limma_res$gene = rownames(limma_res)
    mat_for_limma$gene = rownames(mat_for_limma)
    limma_res = merge(limma_res,mat_for_limma,by='gene')
    
    limma_res = limma_res[order(limma_res$logFC),]
    return(limma_res)


}

do_deseq2 <- function(deseq2_exp_mat, condition_table) {
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
    # output
    resdata$gene = resdata$row.names
    return(resdata)
}

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GOplot)


# input: 
## gene_df:two columns:gene,score_for_rank maybe logFC
# http://yulab-smu.top/clusterProfiler-book/chapter12.html#goplot
# gene_df :symbol and gmt: also symbol
GSEA_enrich <- function(gene_df, gmtfile, pvalueCutoff=0.05){

    c5 <- read.gmt(gmtfile)
    #g_ls = gene_df[order(gene_df$score_for_rank),]$gene
    g_ls = gene_df[order(gene_df$score_for_rank,decreasing = TRUE),]$score_for_rank
    names(g_ls) = gene_df[order(gene_df$score_for_rank,decreasing = TRUE),]$gene

    gsea <- GSEA(g_ls, TERM2GENE = c5, 
                verbose = FALSE, pvalueCutoff = pvalueCutoff,
                pAdjustMethod = 'BH', eps = 1e-10,
                maxGSSize = 500, minGSSize = 5); head(gsea)

    return(gsea)
}

##### dot plot
# plot enrich result(good for gsea result)
## df_for_plot: Description; effect_size; pvalue; set_size
##
## 
## df_for_plot=data.frame(Description=gsea_df$Description, effect_size=gsea_df$NES, pvalue=gsea_df$p.adjust, set_size=gsea_df$setSize)
## text_list=list(pvalue='adj.pval',setSize='setSize',effect_size='NES',y='',title='')
##
## plot_enrich_result(df_for_plot, text_list=text_list)
plot_enrich_result <- function(df_for_plot, text_list=text_list, param_list=NULL, filename=NULL, width=10, height=6) {
    # filter df_for_plot

    #
    df_for_plot$Description = factor(df_for_plot$Description, levels = as.vector(df_for_plot[order(df_for_plot$effect_size),][['Description']]))

    pp = ggplot(df_for_plot,aes(x=effect_size,y=Description))
    #pp + geom_point()
    #pp + geom_point(aes(size=setSize))
    pbubble = pp + geom_point(aes(size=set_size,color=pvalue))
    #pbubble + scale_colour_gradient(low="green",high="red")
    pr = pbubble + scale_colour_gradient(low="green",high="red") 
    pr = pr + labs(color=text_list$pvalue,size=text_list$setSize,x=text_list$effect_size,y=text_list$y,title=text_list$title)
    pr = pr + theme_bw()
    plot(pr)
    
    if (is.null(filename)!=TRUE){
        ggsave(filename, dpi=100, width=width, height=height)
    }

}


# =================

# 表达数据处理
raw_exp_mat = readxl::read_excel('./GSE79888_GEO_value_GCB.DLBCL_RNAseq.xlsx')[c(-1,-2)]
agg_exp_mat = aggregate(.~ gene_name, raw_exp_mat, mean)%>%column_to_rownames('gene_name')

# =================


gmtfile_gobp = '~//my_project/data/pathway/c5.go.bp.v7.2.symbols.gmt'
gmtfile_kegg = '~//my_project/data/pathway/c2.cp.kegg.v7.2.symbols.gmt'


result_ls = list()
# ============================
# ============================
comp = 'Farage'
condition_table = c(rep('Farage.BCL6',3), rep('Farage.DMSO',3))
mat_for_limma = agg_exp_mat[,c(13:18)]



#======= 一下是流程
limma_res = do_limma(mat_for_limma=mat_for_limma, condition_table=condition_table) # 可用deseq2
write_tsv(limma_res, file=sprintf('./%s_limma_res_df.tsv',comp))


## enrichment

gene_df = data.frame(gene=limma_res$gene, score_for_rank=limma_res$t) # 可用deseq2

# KEGG
# gmtfile = '~//my_project/data/pathway/c5.go.bp.v7.2.symbols.gmt'
# gmtfile = '~//my_project/data/pathway/c2.cp.kegg.v7.2.symbols.gmt'

gsea_res_kegg = GSEA_enrich(gene_df, gmtfile_kegg, pvalueCutoff=0.05)
gsea_df = as.data.frame(gsea_res_kegg)


df_for_plot=data.frame(Description=gsea_df$Description, effect_size=gsea_df$NES, pvalue=gsea_df$p.adjust, set_size=gsea_df$setSize)
text_list=list(pvalue='adj.pval',setSize='setSize',effect_size='NES',y='',title='')

plot_enrich_result(df_for_plot, text_list=text_list, param_list=NULL, filename=sprintf('./%s_kegg_gsea.pdf',comp), width=10, height=6)
write_tsv(gsea_df, file=sprintf('./%s_kegg_gsea_df.tsv',comp))



# gobp
#gmtfile = '~//my_project/data/pathway/c5.go.bp.v7.2.symbols.gmt'
#gmtfile = '~//my_project/data/pathway/c2.cp.kegg.v7.2.symbols.gmt'

gsea_res_gobp = GSEA_enrich(gene_df, gmtfile_gobp, pvalueCutoff=0.05)

gsea_df = as.data.frame(gsea_res_gobp)


df_for_plot=data.frame(Description=gsea_df$Description, effect_size=gsea_df$NES, pvalue=gsea_df$p.adjust, set_size=gsea_df$setSize)
text_list=list(pvalue='adj.pval',setSize='setSize',effect_size='NES',y='',title='')
plot_enrich_result(df_for_plot, text_list=text_list, param_list=NULL, filename=sprintf('./%s_gobp_gsea.pdf',comp), width=10, height=6)

write_tsv(gsea_df, file=sprintf('./%s_gobp_gsea_df.tsv',comp))

# ==
sel_gene = limma_res[order(limma_res$adj.P.Val),][c(1:200),]$gene  # 可用deseq2
sel_gene = sapply(strsplit(sel_gene,','),'[[',1)
cat(sel_gene, sep='\n')


result_ls$comp = list(limma_res=limma_res, gsea_res_gobp=gsea_res_gobp, gsea_res_kegg=gsea_res_kegg, sel_gene=sel_gene)  # 可用deseq2













