







dds <- DESeqDataSetFromMatrix(deseq2_exp_mat, 
                              DataFrame(condition_table), 
                              design= ~ condition_table)
dds <- dds[rowSums(counts(dds)) > 1,]
dds2 <- DESeq(dds)
resultsNames(dds2)
# acquire the results using function results(), and assign to res
res <- results(dds2)


