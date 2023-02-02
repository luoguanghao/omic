
if(FALSE){
    # org.Hs.eg.db KEY:
    # ACCNUM,ALIAS,ENSEMBL,ENSEMBLPROT,ENSEMBLTRANS,ENTREZID,ENZYME,EVIDENCE,EVIDENCEALL,GENENAME,GENETYPE,GO,GOALL,IPI,MAP,OMIM,ONTOLOGY,ONTOLOGYALL,PATH,PFAM,PMID,PROSITE,REFSEQ,SYMBOL,UCSCKG,UNIPROT

    library("clusterProfiler")
    library("org.Hs.eg.db")
    library(tidyverse)

    file = '/mnt/d/my_project/zahuo/thg/NASH-转录数据-未标记人类的均为小鼠物种/NASH-转录数据-未标记人类的均为小鼠物种/bulkrnaseq/wholeliver/人类norm_deseq2_whole-gse126848.csv'
    deseq2_df = read_csv(file)

    df2 <- bitr(deseq2_df$Row.names,fromType = "ENSEMBL",toType = c("ENSEMBL","SYMBOL", "GENENAME"),OrgDb = org.Hs.eg.db)
    merge_df = merge(df2,deseq2_df,by.x='ENSEMBL',by.y='Row.names',all.y = TRUE)

    write_tsv(merge_df,file='人类norm_deseq2_whole-gse126848_symbol.csv')

}

