



# DESeq2 #########
deseq2_exp_mat = cbind(exp_mat[,c(16:18)],exp_mat[,c(22:24)])
condition_table = factor(c(rep('dox_off_neg_sgrna_7d',3),rep('dox_on_neg_sgrna_7d',3)))

###
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
resdata = resdata[order(resdata$log2FoldChange),]
# output
head(resdata)
#############





norm_exp_mat = resdata[,c(-1,-2,-3,-4,-5,-6,-7)]
rownames(norm_exp_mat) = resdata$Row.names

upg = resdata$Row.names[(length(resdata$Row.names)-10):length(resdata$Row.names)]
downg = resdata$Row.names[1:10]
# post analysis


## heatmap

heat<-norm_exp_mat[c(downg,upg),]  #对前15个上调基因和下调基因做热图

library(pheatmap)

x <- t(scale(t(heat)))  #事实证明用这个做归一化，效果最好

pheatmap(x,filename='heatmap.png',height=10)





## volcano
library(ggplot2)
library(ggrepel)

library(dplyr)

data<-resdata[c(1:7)]
data$sig[data$pvalue>=0.05 | abs(data$log2FoldChange) < 1] <- "black"
data$sig[data$pvalue<0.05 & data$log2FoldChange >= 1] <- "red"
data$sig[data$pvalue<0.05 & data$log2FoldChange <= -1] <- "green"

input <- data
mark_gene<-c(downg,upg)
volc = ggplot(input, aes(log2FoldChange, -log10(pvalue))) +

    geom_point(aes(col=sig)) + scale_color_manual(values=c("black","green", "red")) +
    geom_point(data=input[input$Row.names %in% mark_gene,], aes(log2FoldChange, -log10(pvalue)), colour="blue", size=2) +
    geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
    geom_vline(xintercept = c(log2(2),-log2(2)),lty=4,lwd=0.6,alpha=0.8)+
    geom_text_repel(data=input[input$Row.names %in% mark_gene,], aes(label=Row.names))

volc




## enrichment

library(clusterProfiler)
library(pathview)

# rt <- read.table("dif.csv",header = T,sep = ",",stringsAsFactors = F)
rt = resdata
colnames(rt)[1] = 'SYMBOL'
GeneSymbol <- rt$SYMBOL
gene.symbol.eg<- id2eg(ids=GeneSymbol, category='SYMBOL', org='Hs',na.rm = F) #如果是小鼠，就是Mm
merg<-merge(gene.symbol.eg,rt,by="SYMBOL")
mer <- subset(merg,merg$ENTREZID != "NA")
geneFC <- mer$log2FoldChange
gene <- mer$ENTREZID
names(geneFC) <- gene

kkd <- enrichGO(gene = gene,"org.Hs.eg.db",ont = "BP",qvalueCutoff = 0.05) #如果是小鼠，就是Mm
write.table(kkd@result,file = "GO.xls",sep = "\t",quote = F,row.names = F)
barplot(kkd,drop=T,showCategory = 12)
cnetplot(kkd,categorySize = "geneNum",foldChange = geneFC)

kk <- enrichKEGG(gene = gene,organism = "human",qvalueCutoff = 0.000001)   #如果是小鼠，就是mouse
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
barplot(kk,drop=T,showCategory = 12)
cnetplot(kk,categorySize = "geneNum",foldChange = geneFC)

## pathview

keggxls<-read.table("KEGG.xls",sep = "\t",header = T)

dir.create("./pathview")
setwd("./pathview") #新建pathview文件夹，存放KEGG关系图

keggxls<-subset(keggxls,keggxls$p.adjust<0.05)
for (i in keggxls$ID){
    pv.out <- pathview(gene.data = geneFC, pathway.id = i, species = "hsa", out.suffix = "pathview") #如果是小鼠mmu，看GSE70410
}




## ppi !!!!
ppi<-read.table("string_interactions0.5.tsv",header = T,sep = "\t",stringsAsFactors = F)
ppi_attr1 <- rt[rt$SYMBOL %in% unique(c(ppi[,1],ppi[,2])),]  # rt即是dif
ppi_attr <- merge(ppi_attr1,as.data.frame(table(c(ppi[,1],ppi[,2]))),by.x="SYMBOL",by.y="Var1")
write.csv(ppi_attr,file = "ppi_attrabute.csv",quote = F,row.names = F)
write.table(ppi_attr,file = "ppi_attrabute.txt",quote = F,sep = "\t",row.names = F)







