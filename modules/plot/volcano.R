

#
# http://202.127.26.99:6424/notebooks/my_home/my_project/other_analysis/zhou/zhou_chk1_depmap.ipynb
# http://202.127.26.99:6424/notebooks/my_home/my_project/other_analysis/yxy_volcano/yxy_volcano.ipynb
#
volc = ggplot(df_for_plot, aes(x=`Effect size`, y=-log10(`P-value`),shape=drug)) +
    geom_point() + scale_color_manual(values=c("#4169E1","grey", "#DC143C")) + # 这个和factor顺序有关
    ggtitle("Wee1")+# geom_point(data=data[data$gene %in% marked_genes,], aes(log2FoldChange, -log10(padj)), colour="blue", size=2) +
    geom_hline(yintercept = -log10(th_p),lty=4,lwd=0.6,alpha=0.8)+
    geom_vline(xintercept = c(th_fc,-th_fc),lty=4,lwd=0.6,alpha=0.8)

volc = volc+geom_text_repel(data=df_for_plot[df_for_plot[['Cancer feature']]=='FLT3_mut',], aes(label=`Cancer feature`), max.overlaps=20) + theme_bw(base_size=15)
# 可以换成 volc = volc+ggrepel::geom_label_repel(data=df_for_plot[df_for_plot[['Cancer feature']]=='FLT3_mut',], aes(label=`Cancer feature`), max.overlaps=20,size=5)
volc = volc + geom_point(data=df_for_plot[df_for_plot[['Cancer feature']]=='FLT3_mut',], aes(x=`Effect size`, y=-log10(`P-value`)), color='red',size=5)
# ggsave('Wee1_plot.pdf', dpi=50, width=8, height=8)


