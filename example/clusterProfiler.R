# R 作图
# https://www.jianshu.com/p/feaefcbdf986
# 使用参看http://n102.yfish.x:8888/notebooks/top/y/Bondi/jupyter/jhw/JHU%E4%BB%A3%E8%B0%A2%E5%B7%AE%E5%BC%82/clusterProfiler.ipynb

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GOplot)

library(dplyr)


# GO ###############################
## bcaa
# data <- read.csv("sites.csv", encoding="UTF-8")
deg_df = read_csv('/y/Bondi/jupyter/jhw/BCAA代谢差异/veh-bcaa-diff_gene.csv')

gene_ls = toupper(deg_df$Row.names) #
file_path = './'
file_name = 'BCAA'
# ===
# input gene_ls:ls of symbol
plot_path = paste(file_path,file_name,'_GO_barplot.pdf',sep='')
table_path = paste(file_path,file_name,'_GOenrich.tsv',sep='')

df2 <- bitr(gene_ls,fromType = "SYMBOL",toType = c("ENTREZID", "GENENAME"),OrgDb = org.Hs.eg.db)
feature_gene_eid = df2[[2]]
kk=enrichGO(gene=feature_gene_eid,OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, ont="ALL", readable=T)

pdf(file=plot_path, width=9, height=7)
bar = barplot(kk, drop=TRUE, showCategory=10, split="ONTOLOGY", color='pvalue') + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
# output: kk: a S4 object, can be tran to data.frame
# write.table(as.data.frame(kk), table_path, row.names=FALSE, col.names=TRUE, sep="\t")
# ===

## jhu
deg_df = read_csv('/y/Bondi/jupyter/jhw/JHU代谢差异/veh-jhu-diff_gene.csv')


gene_ls = toupper(deg_df$Row.names) #

df2 <- bitr(gene_ls,fromType = "SYMBOL",toType = c("ENTREZID", "GENENAME"),OrgDb = org.Hs.eg.db)
feature_gene_eid = df2[[2]]
kk=enrichGO(gene=feature_gene_eid,OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, ont="ALL", readable=T)

## plot graph
pdf(file="./jhu_GO_barplot.pdf", width=9, height=7)
bar = barplot(kk, drop=TRUE, showCategory=10, split="ONTOLOGY", color='pvalue') + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

## jhu
deg_df = read_csv('/y/Bondi/jupyter/jhw/BCAA代谢差异/veh-bcaa-diff_gene.csv')

gene_ls = toupper(deg_df$Row.names) #
df2 <- bitr(gene_ls,fromType = "SYMBOL",toType = c("ENTREZID", "GENENAME"),OrgDb = org.Hs.eg.db)
feature_gene_eid = df2[[2]]
go.BP <- enrichGO(feature_gene_eid, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
plotGOgraph(go.BP)


# KEGG ###############################
## bcaa
deg_df = read_csv('/y/Bondi/jupyter/jhw/BCAA代谢差异/veh-bcaa-diff_gene.csv')

gene_ls = toupper(deg_df$Row.names) #
file_path = './'
file_name = 'BCAA'
# ===
# input gene_ls:ls of symbol
plot_path = paste(file_path,file_name,'_KEGG_barplot.pdf',sep='')
table_path = paste(file_path,file_name,'_KEGGenrich.tsv',sep='')

df2 <- bitr(gene_ls,fromType = "SYMBOL",toType = c("ENTREZID", "GENENAME"),OrgDb = org.Hs.eg.db)
feature_gene_eid = df2[[2]]
kk=enrichKEGG(gene=feature_gene_eid,organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                     minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)

pdf(file=plot_path, width=9, height=7)
bar = barplot(kk, showCategory=30) # there can change
print(bar)
dev.off()
# output: kk: a S4 object, can be tran to data.frame
# write.table(as.data.frame(kk), table_path, row.names=FALSE, col.names=TRUE, sep="\t")

## jhu
deg_df = read_csv('/y/Bondi/jupyter/jhw/JHU代谢差异/veh-jhu-diff_gene.csv')

gene_ls = toupper(deg_df$Row.names) #

df2 <- bitr(gene_ls,fromType = "SYMBOL",toType = c("ENTREZID", "GENENAME"),OrgDb = org.Hs.eg.db)
feature_gene_eid = df2[[2]]
#kk=enrichKEGG(gene=feature_gene_eid,OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, ont="ALL", readable=T)
kk=enrichKEGG(gene=feature_gene_eid,organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                     minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)

pdf(file="./jhu_KEGG_barplot.pdf", width=9, height=7)
dot = dotplot(kk, showCategory=30)
print(dot)
dev.off()

# GSEA ###############################
library(DOSE)
##
deg_df = as.data.frame(read_csv('/y/Bondi/jupyter/jhw/BCAA代谢差异/veh-bcaa-all.csv'))
deg_df$Row.names = toupper(deg_df$Row.names)
rownames(deg_df) = deg_df$Row.names

geneList.sort <- arrange(deg_df, desc(log2FoldChange))$Row.names
eg <- bitr(toupper(geneList.sort), fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")

glist = deg_df[eg$SYMBOL,]$log2FoldChange
names(glist) <- as.character(eg$ENTREZID)

# glist
file_path = './'
file_name = 'BCAA_GOBP'
# ===
# input glist:name is ENTREZID, value is lfc in decreasing
plot_path = paste(file_path,file_name,'_GSEAplot_%s.png',sep='')
table_path = paste(file_path,file_name,'_GSEA.tsv',sep='')

# for general gseGO
gsea_res <- gseGO(glist,OrgDb = org.Hs.eg.db, pvalueCutoff = 0.5)

# for general gsea
# gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
# c5 <- read.gmt(gmtfile)
# gsea <- GSEA(glist, TERM2GENE=c5, verbose=FALSE, pvalueCutoff = 0.8); head(gsea)

for(gene_set_id in 1:10){
    png(file= gsub("%s", gene_set_id, plot_path))
    p = gseaplot2(gsea_res, gene_set_id, title= paste(as.data.frame(gsea_res)[gene_set_id,1],':',as.data.frame(gsea_res)[gene_set_id,2]) )
    print(p)
    dev.off()
    print(gsub("%s", gene_set_id, plot_path))
}

write.table(as.data.frame(gsea_res),table_path,row.names=FALSE,col.names=TRUE,sep="\t")
## =================



deg_df = read_csv('/y/Bondi/jupyter/jhw/BCAA代谢差异/veh-bcaa-all.csv')


deg_df$Row.names = toupper(deg_df$Row.names)
eg <- bitr(toupper(deg_df$Row.names), fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
new_gene_ls = intersect(deg_df$Row.names,eg$SYMBOL)
deg_df = deg_df[which(deg_df$Row.names %in% new_gene_ls),]
deg_df = as.data.frame(deg_df)
rownames(deg_df) = deg_df$Row.names
eg = eg[!duplicated(eg$SYMBOL),]
rownames(eg) = eg$SYMBOL
eg = eg[deg_df$Row.names,]
deg_df$ENTREZID = eg$ENTREZID

geneList <- select(deg_df, ENTREZID, log2FoldChange)
geneList.sort <- arrange(geneList, desc(log2FoldChange))
gene <- geneList.sort$ENTREZID

## feature 1: numeric vector
glist <- geneList[,2]
## feature 2: named vector
names(glist) <- as.character(geneList[,1])
## feature 3: decreasing order
glist <- sort(glist,decreasing = T)

# for general gsea
# gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
# c5 <- read.gmt(gmtfile)
# gsea <- GSEA(glist, TERM2GENE=c5, verbose=FALSE, pvalueCutoff = 0.8); head(gsea)

# for general gseGO
gsea.go <- gseGO(glist,OrgDb = org.Hs.eg.db, pvalueCutoff = 0.5); head(gsea.go)

library(DOSE)
# gseaplot(gsea, 1)
# gseaplot2(gsea, 1)
# gsea.go <- gseGO(glist,OrgDb = org.Hs.eg.db, pvalueCutoff = 0.5); head(gsea.go)

# gene_set_id = 1
for(gene_set_id in 1:10){
    png(file=paste("./GSEA/BCAA_GSEAplot_",gene_set_id,".png",sep=''))
    p = gseaplot2(gsea.go, gene_set_id, title= paste(as.data.frame(gsea.go)[gene_set_id,1],':',as.data.frame(gsea.go)[gene_set_id,2]) )
    print(p)
    dev.off()
    print(paste("./GSEA/BCAA_GSEAplot_",gene_set_id,".png",sep=''))
}
