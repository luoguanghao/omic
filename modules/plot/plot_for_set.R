


## get and show overlap
# some gene_ls combined in list format
# 
# input: all_gene_list, num_overlap
# 获取一个矩阵，每列是一个药，每行是一个基因，描述了每个药都有哪些基因被选中
# output: overlap_gene, maybe 3 overlap or 4 overlap
# 
get_plot_overlap <- function(all_gene_list, num_overlap){

    all_gene = c()
    for(i in names(all_gene_list)){
        all_gene = c(all_gene, all_gene_list[[i]])
    }
    all_gene = unique(all_gene)

    drug_ls = c()
    value_ls = c()
    for(i in names(all_gene_list)){
        drug_ls = c( drug_ls, rep(i,length(all_gene)) )
        value_ls = c( value_ls, all_gene%in%all_gene_list[[i]] )
    }

    df_for_plot = as.data.frame(list( gene=rep(all_gene,4),
                                    drug=drug_ls,
                                    value=value_ls
                                    ))

    df_mean=group_by(df_for_plot%>%filter(value==TRUE), gene) %>% summarise(count = n())
    df_for_plot$gene = factor(df_for_plot$gene, levels = as.vector(df_mean[order(df_mean$count),][['gene']]))

    cols=c('TRUE'='red','FALSE'='black')

    p = df_for_plot%>%ggplot(aes(x=drug,y=gene))+
    geom_tile(aes(fill=value),color="white",size=1)+ #color和size分别指定方块边线的颜色和粗细
    scale_x_discrete("",expand = c(0,0))+ #不显示横纵轴的label文本；画板不延长
    scale_y_discrete("",expand = c(0,0))+
    coord_equal(0.1) +
    scale_fill_manual(values = cols)+ #指定自定义的颜色
    theme(
        axis.text.x.bottom = element_text(size=10),axis.text.y.left = element_text(size = 12), #修改坐标轴文本大小
        axis.ticks = element_blank(), #不显示坐标轴刻度
        legend.title = element_blank() #不显示图例title
    )
    plot(p)
    
    overlap_gene = as.vector( (df_mean%>%filter(count>=num_overlap))$gene )

    return(overlap_gene)

}


































