library(dplyr)
library(DESeq2)
library(tidyverse)

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GOplot)
library(readr)
library(readr)



ENTREZID2SYMBOL <- function(x){
    res = paste(bitr(strsplit(x,'/')[[1]], 
               fromType="ENTREZID", toType=c("SYMBOL","ENSEMBL"), 
               OrgDb="org.Hs.eg.db")$SYMBOL,
         collapse ='//')
    return(res)
    }


# import data
exp_mat = read.csv('./CC-885_exp.csv')

exp_mat <- aggregate(.~ gene_id, exp_mat, mean)
exp_mat = round(exp_mat%>%column_to_rownames('gene_id'))



# prepare for deseq2



# DESeq2 #########
condition_table = factor(c(rep('sen',6),rep('unsen',6)))
deseq2_exp_mat = exp_mat
#
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
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
# output
head(resdata)
#############




# GSEA #################
file_name_prefix = 'kegg'
deg_df = resdata
gmtfile <- '/y/Archive/Bondi/data/pathway/c2.cp.kegg.v7.4.entrez.gmt'
#
deg_df$Row.names = toupper(deg_df$Row.names)
eg <- bitr(toupper(deg_df$Row.names), fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db"); head(eg)

new_gene_ls = intersect(deg_df$Row.names,eg$SYMBOL)
deg_df = deg_df[which(deg_df$Row.names %in% new_gene_ls),]
deg_df = as.data.frame(deg_df)
rownames(deg_df) = deg_df$Row.names
eg = eg[!duplicated(eg$SYMBOL),]
rownames(eg) = eg$SYMBOL
eg = eg[deg_df$Row.names,]
deg_df$ENTREZID = eg$ENTREZID


geneList <- dplyr::select(deg_df, ENTREZID, log2FoldChange); head(geneList) ### can change
geneList.sort <- arrange(geneList, desc(log2FoldChange)); head(geneList.sort) ### can change
gene <- geneList.sort$ENTREZID

c5 <- read.gmt(gmtfile)
# enrich <- enricher(gene, TERM2GENE=c5); head(enrich)
## feature 1: numeric vector
glist <- geneList[,2];head(glist)
## feature 2: named vector
names(glist) <- as.character(geneList[,1]);head(glist)
## feature 3: decreasing order
glist <- sort(glist,decreasing = T); head(glist)

gsea <- GSEA(glist, TERM2GENE=c5, verbose=FALSE, pvalueCutoff = 0.8,minGSSize=5); head(gsea)
## other enrich
## 1 GO ##
#kk <- gseGO(glist,ont = "BP",OrgDb = org.Hs.eg.db)
#gseGO_result<-as.data.frame(kk@result)
#View(gseGO_result)
## 2 KEGG ##
# kk <- gseKEGG(geneList, organism = "mmu",pvalueCutoff = 1)
# gseKEGG_result<-as.data.frame(kk@result)
# gseaplot2(kk, "mmu04060")


##### output
gsea_df = as.data.frame(gsea)
## plot

gsea_df$Description = factor(gsea_df$Description, levels = as.vector(gsea_df[order(gsea_df$enrichmentScore),][['Description']]))

pp = ggplot(gsea_df,aes(x=NES,y=Description))
#pp + geom_point()
#pp + geom_point(aes(size=setSize))
pbubble = pp + geom_point(aes(size=setSize,color=p.adjust))
#pbubble + scale_colour_gradient(low="green",high="red")
pr = pbubble + scale_colour_gradient(low="green",high="red") + labs(color=expression(p.adjust),size="Gene number",x="Enrihment Score",y="Transcription_factor_target_pathway",title="KEGG pathway enrichment")
pr = pr + theme_bw()
ggsave("kegg_dotplot.png",dpi=100,width=10,height=6)

gsea_df$core_enrichment = lapply(gsea_df$core_enrichment,ENTREZID2SYMBOL)

gsea_df = gsea_df%>%filter(p.adjust<0.05) # there should change
write.csv(DataFrame(gsea_df), file= sprintf('./%s_signif_pathway.csv',file_name_prefix))
#################


















