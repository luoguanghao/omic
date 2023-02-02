# 画突变与不突变下药敏
# /d/my_project/jn/5-10分析 仿照AML的分析_fromyfish/FIG4 plot2 boxplot.ipynb
# 



jn_plot4 <- function(res_df,result_dir,filename,save_plot=TRUE){
    ###
    # res_df:(lab_id);ic50;flag;drug;pandel
    # 
    #-i res_df
    res_df[2:4] <- lapply( res_df[2:4], factor)
    res_df$ic50 = res_df$ic50  ###############
    p = ggplot(data = res_df, aes(x = drug, y = ic50, fill = factor(flag)))
    #%R p = p + geom_violin(aes(color=flag))
    p = p + geom_boxplot(outlier.colour = NA,linetype="dashed") ## 设置虚线
    # %R p = p + stat_boxplot(aes(ymin=..lower..,ymax=..upper..,fatten = NULL))  # 去掉须线的纯框
    p = p + geom_boxplot(aes(ymin=..lower.., ymax=..upper..),fatten = NULL, outlier.size = 0.5)
    p = p + stat_summary(fun.y=mean, geom="point",shape=8,position = position_dodge(0.8),color='black')
    # %R p = p + stat_boxplot(geom = "errorbar",aes(x = drug,y=ic50,color = 'black'),width=0.2)
    # %R p = p + stat_boxplot(geom = "errorbar",aes(ymax=..ymin.., fill = factor(flag)),width=0.2)
    # %R p = p + coord_fixed(ratio=1.5)
    p = p + facet_grid(.~pandel, scales = "free",space="free")

    #%R plot(p+labs(title='PI3K/MTOR resistance in CEBPA mutation PDC sample')+theme(text=element_text(size=15,family="serif")))
    # %R p=p+labs(title='PI3K/MTOR resistance in CEBPA mutation PDC sample')
    # %R p=p+theme(text=element_text(size=15,family="serif"))
    p=p+scale_fill_manual("flag",values = mycolors)+scale_color_manual("flag",values = mycolors)
    p = p + theme(text=element_text(size=12,  family="Arial",face = "bold"),
                    axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5), 
                    panel.background = element_blank(),axis.line = element_line(colour="black"),
                    legend.key = element_rect(fill = "white"))
    if(save_plot){
    ggsave(sprintf("%s/%s.png",result_dir,filename), units="in", 
            dpi=100, width=12, height=3.5, device="png")
    }
    return(p)
}







