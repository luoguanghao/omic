library(pheatmap)


## 画heatmap #
# input: df_for_plot, anno_nc_df
## df_for_plot:matrix, context is numeric;colname&rownames
#
## anno_nc_df: rownames is the name in this axis, type col describe the type of each sample
## anno_nc_df = data.frame(type=c(rep('sen',length(sen_sample)),rep('unsen',length(unsen_sample))),row.names=c(sen_sample,unsen_sample))
pheatmap(
    df_for_plot,
    cluster_cols = FALSE,
    #cluster_rows = FALSE,
    scale = 'row',
    color = colorRampPalette(c("green","black", "red"))(50),
    # filename='./cencluster_hmap1.png',
    cellwidth = 8, cellheight = 18,border=FALSE,
    #cutree_cols=2, 
    annotation_col  = anno_nc_df,
    main = drug
    )







