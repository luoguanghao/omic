library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggbeeswarm)



## drug sen data
drug_sen_df = read.csv('/mnt/d/my_project/data/VIZOME/nature_aml_drug_sen.tsv',sep='\t')
drug_sen_df$lab_id = lapply(drug_sen_df$lab_id, function(x) str_replace_all(paste('X',x,sep=''),'-','.'))
drug_sen_df$lab_id = as.factor(as.character(drug_sen_df$lab_id))
drug_ls = drug_sen_df[['inhibitor']][!duplicated(drug_sen_df[['inhibitor']])]

our_drugs = c('GDC-0941','BEZ235','INK-128','MK-2206','PI-103')


## mut data
maf_df = read.csv('/mnt/d/my_project/data/VIZOME/LAML_nature_PDC.maf',sep='\t')
gene_count_df = as.data.frame(table(maf_df$Hugo_Symbol))
gene_count_df = gene_count_df[order(gene_count_df$Freq,decreasing =TRUE),]

top_gene = as.vector(gene_count_df$Var1[1:20])

mut_sample_ls = as.vector( (maf_df %>% filter(Hugo_Symbol==tmp_mut))['Tumor_Sample_Barcode'] )[['Tumor_Sample_Barcode']] ### <<<<
mut_sample_ls = gsub('-','.',paste('X',mut_sample_ls,sep=''))


drug = our_drugs[1] #### <<<

jn_plot <- function(drug, mut_sample_ls, top_gene, maf_df, drug_sen_df){
    df_for_plot = 'NA'

    for(i in 1:length(top_gene)){
        tmp_mut = top_gene[i]

        mut_sample_ls = as.vector( (maf_df %>% filter(Hugo_Symbol==tmp_mut))['Tumor_Sample_Barcode'] )[['Tumor_Sample_Barcode']] ### <<<<
        mut_sample_ls = gsub('-','.',paste('X',mut_sample_ls,sep=''))
        
        sp_drug_sen_df = drug_sen_df%>%filter(inhibitor==drug)
        tmp_df_for_plot = sp_drug_sen_df%>%filter(lab_id%in%mut_sample_ls)
        tmp_df_for_plot$mut = tmp_mut
        if(i==1){
            df_for_plot = tmp_df_for_plot
        }else{
            df_for_plot = rbind(df_for_plot,tmp_df_for_plot)
        }
    }
    
    # get bar_df
    attach(df_for_plot)
    mean_df<-aggregate(df_for_plot[c('mut','ic50','auc')],by=list(mut),FUN=mean,na.rm=TRUE)
    mean_df$mut = mean_df$Group.1

    mean_df$mut = factor(mean_df$mut,level=as.vector(mean_df$mut[order(mean_df$auc)]))
    df_for_plot$mut = factor(df_for_plot$mut,level=as.vector(mean_df$mut[order(mean_df$auc)]))

    # plot
    p = ggplot()
    p = p + geom_bar(data=mean_df,aes(x=mut,y=auc),stat = "identity",fill='DarkBlue',width=0.7)
    p = p + geom_beeswarm(data=df_for_plot,aes(x=mut, y=auc),color='red',size=1,dodge.width=0.8) 
    p = p+xlab('')
    p = p+ylab('AUC')
    p = p + scale_y_continuous(expand=c(0,0))
    p = p + theme(text=element_text(size=16, face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5),
                panel.background=element_blank(),,axis.line = element_line(colour="black"))
    ggsave(sprintf('fig7_%s.pdf',drug), dpi=50, width=12, height=6)

    return(p)

}


signif_test_for_mut <- function(drug, mut_sample_ls, top_gene, maf_df, drug_sen_df){





}






