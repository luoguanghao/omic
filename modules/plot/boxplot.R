# ggplot and then ggsave()
# 
# ggsave(): 2condition boxplot just width=5,height=5
#
#
#
#
library(maftools)
library(ggpubr)
library(ggplot2)
# library(ploty)



################
# 分组barplot
################
# input: a long frame: sample;drugsen;mutwt
# process: calculate mean,sd,count,s.e.m
#          calculate test for every drug
# output: boxplot, test result

get_bar_df <- function(df_for_plot){

    attach(df_for_plot)
    bar_df<-aggregate(df_for_plot[c('mut','ic50','auc')],by=list(mut),FUN=mean,na.rm=TRUE)
    bar_df$mut = bar_df$Group.1

    sd_df<-aggregate(df_for_plot[c('mut','ic50','auc')],by=list(mut),FUN=sd,na.rm=TRUE)
    sd_df$mut = sd_df$Group.1
    count_df<-as.data.frame(table(df_for_plot$mut))
    colnames(count_df) = c('gene','count')
    
    bar_df$ic50_sem = sd_df$ic50/(count_df$count)**0.5
    bar_df$auc_sem = sd_df$auc/(count_df$count)**0.5
    
    
    bar_df$mut = factor(bar_df$mut,level=as.vector(bar_df$mut[order(bar_df$auc)]))
    # sd_df$mut = factor(sd_df$mut,level=as.vector(bar_df$mut[order(bar_df$auc)]))
    df_for_plot$mut = factor(df_for_plot$mut,level=as.vector(bar_df$mut[order(bar_df$auc)]))

    
    bar_df$ic50_sd = sd_df$ic50
    bar_df$auc_sd = sd_df$auc

    return(bar_df)

}

my_barplot <- function(df_for_plot, label_list, para_list, out_name) { # <<<


    bar_df <- get_bar_df(df_for_plot)
    #mean_df



    p = ggplot()
    p = p + geom_errorbar(data=bar_df,aes(x=mut,ymin = auc - auc_sem, ymax = auc + auc_sem), width = 0.2)
    p = p + geom_bar(data=bar_df,aes(x=mut,y=auc),stat = "identity",fill='DarkBlue',width=0.7)
    p = p + geom_beeswarm(data=df_for_plot,aes(x=mut, y=auc),color='red',size=1,dodge.width=0.8) 
    # title ....
    p = p+xlab('')
    p = p+ylab('AUC')
    #p = p+ylim(0, 310)
    p = p + scale_y_continuous(limit=c(0,310),expand = c(0,0))
    p = p + theme(text=element_text(size=16, face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5),
                panel.background=element_blank(),axis.line = element_line(colour="black"))
    ggsave(sprintf('fig7_%s.pdf',drug), dpi=50, width=12, height=6) # <<<<
    plot(p)

    return(p)

}



##########################################
# 融合方法，融合vol图，boxplot，jitter图  USEFUL!!!
##########################################
# df_for_plot: x_val,y_val,(c_val),other_info
# label_list: $ylab/xlab/labs  list(ylab=, labs=)
# para_list
# out_name
my_boxplot <- function(df_for_plot, x='x_val', y='y_val', label_list=list(labs='plot',ylab='',xlab=''), para_list=NA, out_name=NA) { # <<<
    mycolors = c('#FF0000','#225EA8')
    df_for_plot[[x]] = as.factor(df_for_plot[[x]])
    
    p = ggplot(data=df_for_plot,aes(x=get(x), y=get(y)))
    # p = p + geom_violin(aes(fill=x_val),color='black') # fill maybe c_val
    p = p + geom_boxplot(aes(fill=get(x)),color='black')

    p = p + geom_jitter(fill='white',width =0.1,shape = 21,size=1.5)
    # p = p + geom_beeswarm(data=df_for_plot,aes(x=mut, y=auc),color='red',size=1,dodge.width=0.8) 

    # p = p+stat_summary(fun.y=mean, geom="point")
    p = p + stat_compare_means()
    # title ....
    p = p+labs(title=label_list$labs)
    p = p+ylab(label_list$ylab) + xlab(label_list$xlab)

    p=p+scale_fill_manual(x,values = mycolors) + scale_color_manual(y,values = mycolors) # fill maybe c_val

    p=p + theme(text=element_text(size=20,  family="Arial"),
            axis.text.y=element_text(size=20),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=20),
            panel.background = element_blank(),
            axis.line = element_line(colour="black"),
            legend.key = element_rect(fill = "white"))

    if(is.na(out_name)==FALSE){
        ggsave(sprintf(out_name), dpi=50, width=12, height=6)
    } # <<<<<<<<
    plot(p)

    return(p)

}




old_my_boxplot <- function(df_for_plot, label_list, para_list=NA, out_name=NA) { # <<<
    mycolors = c('#FF0000','#225EA8')

    p = ggplot(data=df_for_plot,aes(x=x_val, y=y_val))
    # p = p + geom_violin(aes(fill=x_val),color='black') # fill maybe c_val
    p = p + geom_boxplot(aes(fill=x_val),color='black')

    p = p + geom_jitter(fill='white',width =0.1,shape = 21,size=1.5)
    # p = p + geom_beeswarm(data=df_for_plot,aes(x=mut, y=auc),color='red',size=1,dodge.width=0.8) 

    # p = p+stat_summary(fun.y=mean, geom="point")
    p = p + stat_compare_means()
    # title ....
    p = p+labs(title=label_list$labs)
    p = p+ylab(label_list$ylab)

    p=p+scale_fill_manual('x_val',values = mycolors) + scale_color_manual('x_val',values = mycolors) # fill maybe c_val

    p=p + theme(text=element_text(size=20,  family="Arial",face='bold'),
            axis.text.y=element_text(size=20),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,size=20),
            panel.background = element_blank(),axis.line = element_line(colour="black"),legend.key = element_rect(fill = "white"))

    if(is.na(out_name)==FALSE){ggsave(sprintf(out_name), dpi=50, width=12, height=6)} # <<<<<<<<
    plot(p)

    return(p)

}

if(FALSE){# 区分形状的violin plot：输入data.frame: cell_line;concentration;CCC_124;solid_blood
    ggpubr::ggviolin(comb_drug_df,x='concentration',y='CCC_124',color='concentration',add='boxplot',size=0.2)+ggpubr::stat_compare_means()+
    geom_jitter(data=comb_drug_df,aes(x=concentration,y=CCC_124,shape=solid_blood,color=concentration),size=3,width=0.2)
    ggpubr::compare_means(CCC_124~concentration,comb_drug_df,p.adjust.method = 'fdr')
}

# 分panel的boxplot
# http://localhost:8888/notebooks/d/my_project/hl/changhai%E9%AA%8C%E8%AF%81/new_8_NPM1_blast_gene_ch_newsample.ipynb
if(0){

    p <- ggboxplot(df_for_plot, 
              x = "npm1_stat", y = "ic50",
              color = "npm1_stat", palette = "nature",
              add = "jitter") + 
    stat_compare_means() +
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) +
    facet_wrap(~drug, strip.position="bottom")
p
}
















