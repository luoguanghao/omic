library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("GOplot")
library(ggthemes)


## over-representation #
#
# type='enrich' for self design enrich, input is the gene_symbol
# http://yulab-smu.top/clusterProfiler-book/chapter12.html#goplot
if(FALSE){
#
# keytypes() 看支持转换类型
#
#
#PATH_DATA = '~/my_project/data/'
gmtfile = list(
   tft=paste('/home/lgh/my_project/data/pathway/c3.tft.v7.5.1.symbols.add_cebpa.gmt',sep=''),
   kegg=paste('/home/lgh/my_project/data/pathway/c2.cp.kegg.v7.5.1.symbols.gmt',sep=''),
   reactome=paste('/home/lgh/my_project/data/pathway/c2.cp.reactome.v7.5.1.symbols.gmt',sep='')
)
}
#
do_ora_enrich <- function(gene_ls, gmt_file, pvalueCutoff=0.05,save_dir=NULL, project_name=NULL){
    res_list = list()
    print(length(gene_ls))
    #cat(sprintf(''))
    for(gmt in names(gmt_file)){
        c5 <- read.gmt(gmt_file[[gmt]])
        kk <- enricher(gene_ls, TERM2GENE=c5, 
                    pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff, 
                    minGSSize = 5, maxGSSize = 500)

        # kk = as.data.frame(kk)
        res_list[[gmt]] = kk
        if(!is.null(save_dir)){
            write_tsv(data.frame(kk),
                    file=sprintf('%s/%s_ora_%s.tsv',save_dir,project_name,gmt))
        }
        # print(sprintf('## enriched %s %s items',dim(data.frame(kk)%>%filter(p.adjust<0.05))[1],gmt ))
    }
    
    return(res_list)
}

OR_enrich <- function(gene_ls, type='enrichGO', gmtfile=NULL, pvalueCutoff=0.05,ont='ALL'){
    
    if(type=='enrichGO'){
        df2 <- bitr(gene_ls,fromType = "SYMBOL",toType = c("ENTREZID", "GENENAME"),OrgDb = org.Hs.eg.db)
        feature_gene_eid = df2[[2]]
        kk=enrichGO(gene=feature_gene_eid,OrgDb=org.Hs.eg.db, pvalueCutoff=pvalueCutoff, ont=ont, readable=T)

        # return(as.data.frame(kk))
    }
    if(type=='enrichKEGG'){
        df2 <- bitr(gene_ls,fromType = "SYMBOL",toType = c("ENTREZID", "GENENAME"),OrgDb = org.Hs.eg.db)
        feature_gene_eid = df2[[2]]
        kk <- enrichKEGG(gene=feature_gene_eid, organism='hsa', pvalueCutoff=pvalueCutoff)

        # return(as.data.frame(kk))
    }
    if(type=='enrich'){

        # gmtfile <- '/mnt/d/my_project/data/pathway/c2.cp.kegg.v7.2.symbols.gmt'
        c5 <- read.gmt(gmtfile)

        kk <- enricher(gene_ls, TERM2GENE=c5, 
                    pAdjustMethod = "BH", pvalueCutoff = pvalueCutoff, 
                    minGSSize = 5, maxGSSize = 500)

        # return(as.data.frame(kk))
    }
    kk = as.data.frame(kk)
    return(kk)

}




# input: 
## gene_df:two columns:gene,score_for_rank maybe logFC
# http://yulab-smu.top/clusterProfiler-book/chapter12.html#goplot
# gene_df :symbol and gmt: also symbol
# type: path , df
GSEA_enrich <- function(gene_df, gmtfile, pvalueCutoff=0.05, type='path'){

    tmp_gene_df = na.omit(gene_df)
    if(dim(tmp_gene_df)[1]!=dim(gene_df)[1]){
        message(sprintf('remove na row, from %s to %s',dim(gene_df)[1],dim(tmp_gene_df)[1]))
    }
    gene_df = tmp_gene_df

    colnames(gene_df) = c('gene','score_for_rank')
    if(type=='path')
        c5 <- read.gmt(gmtfile)
    else
        c5 <- gmtfile

    #g_ls = gene_df[order(gene_df$score_for_rank),]$gene
    g_ls = gene_df[order(gene_df$score_for_rank,decreasing = TRUE),]$score_for_rank
    names(g_ls) = gene_df[order(gene_df$score_for_rank,decreasing = TRUE),]$gene

    gsea <- GSEA(g_ls, TERM2GENE = c5, 
                verbose = FALSE, pvalueCutoff = pvalueCutoff,
                pAdjustMethod = 'BH', eps = 1e-10,
                maxGSSize = 500, minGSSize = 1); head(gsea)

    return(gsea)
}

# do gesGO
if(FALSE){


    g_ls = bcaa_diff_df[order(bcaa_diff_df$stat,decreasing = TRUE),]$stat
    names(g_ls) = toupper(bcaa_diff_df[order(bcaa_diff_df$stat,decreasing = TRUE),]$Row.names)
    df2 <- bitr(names(g_ls),fromType = "SYMBOL",toType = c("ENTREZID", "GENENAME"),OrgDb = org.Hs.eg.db)
    names(g_ls) = df2[[2]]


    bcaa_gsego_res = gseGO(g_ls,ont='ALL',OrgDb=org.Hs.eg.db,maxGSSize = 500, minGSSize = 2)
    bcaa_gsego_df = as.data.frame(bcaa_gsego_res);head(bcaa_gsego_df)


}

# make expmat into GSEA format
# what is that format
if(FALSE){

    
}


#############################
# KSEA
#############################
# pho_df
# stat_res: 'gene','p',effect_name
# effect_name: 'r' or 'fc_median'
# t: sample*gene -> gene*sample
do_KSEA <- function(pho_df,stat_res,effect_name,p_name='p',t=FALSE,unlog=FALSE){
    library(KSEAapp)
    if(t){
        pho_df_t = data.frame(t(pho_df[,-1]))
        colnames(pho_df_t) = pho_df[[1]]
        pho_df = pho_df_t%>%rownames_to_column('site')
    }        
    # return(pho_df)
    pre_px = data.frame(
            key = pho_df[[1]],
            Protein = rep('NULL',nrow(pho_df)),
            Gene=sapply(strsplit(pho_df[[1]],'_'),'[[',1),
            Peptide = rep('NULL',nrow(pho_df)),
            Residue.Both=sapply(strsplit(pho_df[[1]],'_'),'[[',2)
        )

    my_px = merge(pre_px,stat_res[,c('gene',p_name,effect_name)],by.x='key',by.y='gene')[,-1]
    colnames(my_px)[5:6] = c('p','FC')
    if(unlog){my_px$FC = 2**my_px$FC}
    #my_px$FC = 2**(my_px$FC)
    #return(list(pre_px,stat_res,my_px))
    

    ## 处理Inf值0值na值问题

    my_px = na.omit(my_px)

    if(sum(my_px$FC==Inf)>0){
        my_px[my_px$FC==Inf,]$FC = max(
                    my_px[my_px$FC!=Inf,]$FC)
    }

    if(sum(my_px$FC==0)>0){
        my_px[my_px$FC==0,]$FC = min(
                    my_px[my_px$FC!=0,]$FC)
    }

    scores = KSEA.Scores(KSData, my_px, NetworKIN=FALSE)
    KSEA.Barplot(KSData, my_px, NetworKIN=FALSE, m.cutoff=2,p.cutoff=0.1,export=FALSE)

    scores = scores[order(scores$p.value),]

    return(list(ksea_scores=scores,px=my_px))

}


if(FALSE){

    # http://202.127.26.99:6424/notebooks/my_home/my_project/hl/BCL2_%E5%88%86%E6%9E%90/ksea%E5%88%86%E6%9E%90.ipynb

    pre_px = data.frame(
            key = pho_df[[1]],
            Protein = rep('NULL',nrow(pho_df)),
            Gene=sapply(strsplit(pho_df[[1]],'_'),'[[',1),
            Peptide = rep('NULL',nrow(pho_df)),
            Residue.Both=sapply(strsplit(pho_df[[1]],'_'),'[[',2)
        )

    my_px = merge(pre_px,all_res_list$ABT_737[,c('gene','p','r')],by.x='key',by.y='gene')[,-1]
    colnames(my_px)[5:6] = c('p','FC')
    my_px$FC = 2**(my_px$FC)

    scores = KSEA.Scores(KSData, my_px, NetworKIN=FALSE)
    KSEA.Barplot(KSData, my_px, NetworKIN=FALSE, m.cutoff=2,p.cutoff=0.1,export=FALSE)


}







###############################################
# simple plot
## dot plot
## bagua plot
## barplot

## gsea plot
###############################################

##### dot plot
# plot enrich result(good for gsea result)
## df_for_plot: Description; effect_size; pvalue; set_size
##
## gsea_df = as.data.frame(gsea_res)
## df_for_plot=data.frame(Description=gsea_df$Description, effect_size=gsea_df$NES, pvalue=gsea_df$p.adjust, set_size=gsea_df$setSize)
## text_list=list(pvalue='adj.pval',setSize='setSize',effect_size='NES',y='',title='')
##
## plot_enrich_result(df_for_plot, text_list=text_list)
plot_enrich_result <- function(df_for_plot, text_list=text_list, param_list=NULL, filename=NULL, width=10, height=6, show_num=14){
    df_for_plot$pvalue = df_for_plot[[text_list$pvalue]]
    df_for_plot$setSize = df_for_plot[[text_list$setSize]]
    df_for_plot$effect_size = df_for_plot[[text_list$effect_size]]

    df_for_plot = df_for_plot[order(df_for_plot$pvalue),][c(1:show_num),]
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

# func inside the clusterProfiler
#  clusterProfiler::dotplot(kk, x="count", showCategory=20, colorBy="qvalue")






##### bagua plot
# circle_dat(terms, genes)
## terms: A data frame with columns for 'category', 'ID', 'Term', adjusted p-value ('adj_pval') and 'Genes'(sep by ', ')
### df_for_goplot_ls$terms=data.frame(Category=,ID=,Term=,adj_pval=,Genes=)
##
## genes: A data frame with columns for 'ID', 'logFC' (gene id)
### df_for_goplot_ls$genes=data.frame(ID=,logFC=)
# https://wencke.github.io/
#
#
# 八卦图的逻辑是啥?
# z-score:, 排序, 内环高低 
#
if(FALSE){
    # prepare gor gocircle
    ## input: gmt_file; gsea_df:as.data.frame(GSEA_enrich_res); deseq2_df
    gmt_df = read.gmt(gmt_file)

    gmt_list = list()
    for(item in unique(gmt_df$term)){
        tm_gene = intersect(deseq2_df$Row.names,as.vector((gmt_df%>%filter(term==item))$gene))
        gmt_list[[item]] = base::paste(tm_gene,collapse=', ')
    }
    ont_genes_df = data.frame(Term=names(gmt_list),Genes=as.matrix(gmt_list))
    colnames(ont_genes_df) = c('Term','Genes')

    df_for_goplot_ls = list()
    df_for_goplot_ls$terms=data.frame(Category='BP',ID=c(1:length(gsea_df$Description)),Term=gsea_df$Description,adj_pval=gsea_df$p.adjust)
    df_for_goplot_ls$genes=data.frame(ID=deseq2_df$Row.names,logFC=deseq2_df$log2FoldChange)
    df_for_goplot_ls$terms = merge(df_for_goplot_ls$terms,ont_genes_df,by='Term')
    circ <- circle_dat(df_for_goplot_ls$terms, df_for_goplot_ls$genes)
}

if(FALSE){
    circ <- circle_dat(df_for_goplot_ls$terms, df_for_goplot_ls$genes)

    GOBar(subset(circ, category == 'BP'))

    GOBubble(circ, labels = 3)

    pdf('bagua.pdf',width=15,height=10)
    GOCircle(circ)
    dev.off()
}


######## GSEA plot
# input: result of GSEA_enrich, S4 object
if(FALSE){
    id = 'GO:0044241'
    png(sprintf('./%s_gsea_plot.png',str_replace(id,':','_')))
    gseaplot2(bcaa_gsego_res,geneSetID=id, title=name)
    dev.off()

}



##############################################
# workflow
##############################################



#######################
# gmt file
#######################


gmt2list <- function(gmt_file){
    gmt_df = read.gmt(gmt_file)
    gmt_list = list()
    for(item in unique(gmt_df$term)){
        gmt_list[[item]] = as.vector((gmt_df%>%filter(term==item))$gene)
    }
    return(gmt_list)
}













