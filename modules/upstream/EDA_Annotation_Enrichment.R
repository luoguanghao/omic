
#
#
#
library(tidyverse)
library(Rsubread)
library(DESeq2)


SraRunTable_path = '~/raw_data/PRJNA726826/SraRunTable.txt'
exp_mat_path = '../do_subread/PRJNA726826_count_exp_mat.tsv'

# get sample info

SraRunTable = read.table(SraRunTable_path, sep=',',header=1)

anno_table = SraRunTable[,c('Run','Sample.Name')]
colnames(anno_table) = c('Run','SampleName')

anno_table$'group' = sapply(strsplit(anno_table$'SampleName', "\\-"),function(x) paste(x[1:3],collapse='-'))

# get express matrix

exp_mat = read.table(exp_mat_path)

## change colnames

colnames(exp_mat) = sapply(strsplit(colnames(exp_mat), "\\."),'[[',1)
colnames(exp_mat) = c('gene_name', (anno_table%>%filter(Run%in%colnames(exp_mat)[2:10]))$SampleName)



res_ls = list()

## 提取需要的样本，0h-12h
deseq2_exp_mat = (exp_mat%>%column_to_rownames('gene_name'))[( anno_table%>%filter(group%in%c('L02-PAOA-0h','L02-PAOA-12h')) )$SampleName] # << 选取样本
condition_table = as.factor(( anno_table%>%filter(group%in%c('L02-PAOA-0h','L02-PAOA-12h')) )$'group')  # << 选取样本


do_deseq2 <- function(deseq2_exp_mat, condition_table) {
    # ===
    # ii                                       
    dds <- DESeqDataSetFromMatrix(deseq2_exp_mat, 
                                DataFrame(condition_table), 
                                design= ~ condition_table)
    dds <- dds[rowSums(counts(dds)) > 1,]
    dds2 <- DESeq(dds)
    resultsNames(dds2)
    # acquire the results using function results(), and assign to res
    res <- results(dds2)
    # summary(res)
    DEG_exp_df <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
    return(DEG_exp_df)
}

DEG_exp_df <- do_deseq2(deseq2_exp_mat, condition_table)

res_ls$X_0_12 = DEG_exp_df




















