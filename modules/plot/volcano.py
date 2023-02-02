
######
# 1
######
%R res_df$mean_wt = as.numeric(res_df$mean_wt)
%R res_df$mean_mut = as.numeric(res_df$mean_mut)
%R res_df$pvalue = as.numeric(res_df$pvalue)
%R res_df$adjpvalue = as.numeric(res_df$adjpvalue)
%R res_df$FC = as.numeric(res_df$FC)

%R plot_df <- res_df
%R plot_df$label = ''
%R plot_df$label[plot_df$adjpvalue<0.05] = plot_df$lineage_subtype[plot_df$adjpvalue<0.05]

%R p <-ggplot(data=subset(plot_df,adjpvalue!=888), aes(x=log2(FC), y=-log10(adjpvalue)) )
%R p = p + geom_point(alpha=0.8, size=1)  
%R p = p + xlab("log2 fold change") + ylab("-log10 padj")
%R p = p + geom_hline(yintercept=-log10(0.05),linetype=4)
# %R p = p + geom_vline(xintercept=c(-1,1),linetype=4)
%R p = p + geom_text(aes(label = label), size = 3,vjust=-0.5, alpha=0.8,family="Arial",color='red')
%R p = p + theme(text=element_text(size=12,  family="Arial",face = "bold"))
%R plot(p)
#%R ggsave("./volcano.png")

'''add condition
p <-ggplot(data=mydata, aes(x=log2FoldChange, y=-log10(padj), colour=Condition)) + 
    geom_point(alpha=0.8, size=1)  +  xlab("log2 fold change") + ylab("-log10 padj")+
    geom_hline(yintercept=-log10(0.05),linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+scale_color_manual(values=c('up'='red','down'='green','normal'='gray'))
p
'''















