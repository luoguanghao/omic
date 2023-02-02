library(ggplot2)
library(tidyverse)
library(ggpubr)









#
# maf_df
# gene
# our_drugs
# drug_sen_df colnames=c('Tumor_Sample_Barcode'，'inhibitor','auc','ic50')
# result_dir
# aucic50
#
# OUTPUT:list(df_for_plot=df_for_plot, p=p, compare_means_df=compare_means_df)x
jn_plot5 <- function(maf_df,gene,our_drugs,drug_sen_df,aucic50,result_dir,save_plot=FALSE){

    mut_sample_ls = as.vector(
            (maf_df %>% filter(Hugo_Symbol==gene))['Tumor_Sample_Barcode'] )[['Tumor_Sample_Barcode']]
    all_sample_ls = unique(maf_df$Tumor_Sample_Barcode)
    wt_sample_ls = setdiff(all_sample_ls,mut_sample_ls)

    df_for_plot = drug_sen_df%>%filter(inhibitor%in%our_drugs)

    df_for_plot$mut = 'NA'

    df_for_plot[df_for_plot$Tumor_Sample_Barcode%in%mut_sample_ls,]$mut = 'mut'
    df_for_plot[df_for_plot$Tumor_Sample_Barcode%in%wt_sample_ls,]$mut = 'wt'

    df_for_plot = df_for_plot%>%filter(mut!='NA')

    mycolors = c('grey','red')

    p = ggplot(df_for_plot,aes(x = inhibitor, y = aucic50))
    p = p + geom_jitter(aes(fill=mut),color='white',position = position_jitterdodge(jitter.width=0.15),shape=21)
    p = p + scale_fill_manual('mut',values = mycolors)

    p = p + stat_summary(aes(group=mut),fun=mean, geom="point",shape='-',size=10, position = position_dodge(0.8))

    p = p + ylab(aucic50) + xlab('Inhibitors')

    p = p + theme(text=element_text(size=18,face='bold'),axis.text.x = element_text(angle = 40, hjust = 0.5,vjust=0.5),
                axis.title.x = element_text(vjust = -0.5),
                panel.background=element_blank(),axis.line = element_line(colour="black"))
    # print(p)
    if(save_plot){
        ggsave(sprintf('%s/fig5.pdf',result_dir), dpi=50, width=7, height=7)
    }

    compare_means_df = compare_means(formula(sprintf('%s~mut',aucic50)),df_for_plot,group.by='inhibitor',method = "wilcox.test",p.adjust.method='fdr')

    return(list(df_for_plot=df_for_plot, p=p, compare_means_df=compare_means_df))
}



if(FALSE){
## drug sen data
drug_sen_df = read.csv('/mnt/d/my_project/data/VIZOME/nature_aml_drug_sen.tsv',sep='\t')
drug_sen_df$lab_id = lapply(drug_sen_df$lab_id, function(x) str_replace_all(paste('X',x,sep=''),'-','.'))
drug_sen_df$lab_id = as.factor(as.character(drug_sen_df$lab_id))
drug_ls = drug_sen_df[['inhibitor']][!duplicated(drug_sen_df[['inhibitor']])]

our_drugs = c('GDC-0941','BEZ235','INK-128','MK-2206','PI-103')

## mut data
maf_df = read.csv('/mnt/d/my_project/data/VIZOME/LAML_nature_PDC.maf',sep='\t')
maf_df$Tumor_Sample_Barcode = gsub('-','.',paste('X',maf_df$Tumor_Sample_Barcode,sep=''))

jn_plot5(maf_df,'CEBPA',our_drugs,drug_sen_df,'ic50','./result',save_plot=FALSE)

}





