library(ggplot2)
library(tidyverse)
library(ggpubr)



# expression data #################
## from VIZOME
exp_mat = read_tsv('/mnt/d/my_project/data/VIZOME/nature_aml_log2_fpkm.txt')
exp_mat = aggregate(.~ Symbol, exp_mat, mean) #处理重复基因名问题
rownames(exp_mat) = exp_mat$Symbol
exp_mat = exp_mat[-1]




colnames(exp_mat) = paste('X',colnames(exp_mat),sep='')

colnames(exp_mat) = gsub("-",".",colnames(exp_mat) )

## drug sen data
drug_sen_df = read.csv('/mnt/d/my_project/data/VIZOME/nature_aml_drug_sen.tsv',sep='\t')
drug_sen_df$lab_id = lapply(drug_sen_df$lab_id, function(x) str_replace_all(paste('X',x,sep=''),'-','.'))
drug_sen_df$lab_id = as.factor(as.character(drug_sen_df$lab_id))
drug_ls = drug_sen_df[['inhibitor']][!duplicated(drug_sen_df[['inhibitor']])]





# our_drugs = c('Midostaurin','Quizartinib (AC220)','Crenolanib','Gilteritinib (ASP-2215)')
our_drugs = c('GDC-0941','BEZ235','INK-128','MK-2206','PI-103')

#############################################
#############################################
drug = our_drugs[1]

## 看看药敏的分布
ggplot(data=drug_sen_df[drug_sen_df['inhibitor']==drug,]) + 
  geom_histogram(aes(x=auc),bins=30)

df_dot_plot = drug_sen_df[drug_sen_df['inhibitor']==drug,][order(drug_sen_df[drug_sen_df['inhibitor']==drug,]$auc),]
df_dot_plot$lab_id = factor(df_dot_plot$lab_id,level=as.vector(df_dot_plot$lab_id))

plot(x=1:423,y=df_dot_plot$ic50)
ggplot(data=df_dot_plot) + geom_point(aes(x=lab_id, y=auc)) + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))



# analysis
sp_drug_sen_df = drug_sen_df%>%filter(inhibitor==drug)

## divide
### get sen unsen
top_bottom_number = round(dim(sp_drug_sen_df)[1]*0.3)
sp_drug_sen_df = drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],][order(drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],]$auc),]
sp_drug_sen_df$lab_id = factor(sp_drug_sen_df$lab_id,level=as.vector(sp_drug_sen_df$lab_id))
rownames(sp_drug_sen_df) = sp_drug_sen_df$lab_id
sp_sample = intersect(sp_drug_sen_df$lab_id,colnames(exp_mat))

#sen_sample = (sp_drug_sen_df[sp_drug_sen_df$lab_id %in% sp_sample,] %>% filter(auc<182.3626))$lab_id ### should change
#unsen_sample = (sp_drug_sen_df[sp_drug_sen_df$lab_id %in% sp_sample,] %>% filter(auc>232.7781))$lab_id ### should change
sen_cutoff = sp_drug_sen_df[top_bottom_number,]$auc
unsen_cutoff = sp_drug_sen_df[
    dim(sp_drug_sen_df)[1]-top_bottom_number,]$auc


sen_sample = as.vector((sp_drug_sen_df[sp_drug_sen_df$lab_id %in% sp_sample,] %>% filter(auc<sen_cutoff))$lab_id )
unsen_sample = as.vector((sp_drug_sen_df[sp_drug_sen_df$lab_id %in% sp_sample,] %>% filter(auc>unsen_cutoff))$lab_id )



df_fpr_plot = data.frame(sample=colnames(exp_mat),exp=t(exp_mat)[,'CEBPA'],drug_sen=rep('NA',dim(exp_mat)[2]))

df_fpr_plot[sen_sample,]$drug_sen = 'sen'

df_fpr_plot[unsen_sample,]$drug_sen = 'unsen'

df_fpr_plot = df_fpr_plot%>%filter(drug_sen!='NA')



df_fpr_plot$drug_sen = as.factor(df_fpr_plot$drug_sen)

# p = ggplot(df_fpr_plot,aes(x = drug_sen, y = exp,fill=drug_sen)) + geom_boxplot()
p <- ggboxplot(df_fpr_plot, x = 'drug_sen', y = 'exp',add = "jitter",color='drug_sen', font.label=) + stat_compare_means()
# p = p + theme(text=element_text(size=15,  family="Times New Roman"))
plot(p)



#################################
#################################







jn_boxplot <- function(drug,exp_mat,drug_sen_df){
    
    # analysis
    sp_drug_sen_df = drug_sen_df%>%filter(inhibitor==drug)

    ## divide
    ### get sen unsen
    top_bottom_number = round(dim(sp_drug_sen_df)[1]*0.3)
    sp_drug_sen_df = drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],][order(drug_sen_df[drug_sen_df['inhibitor']==our_drugs[1],]$auc),]
    sp_drug_sen_df$lab_id = factor(sp_drug_sen_df$lab_id,level=as.vector(sp_drug_sen_df$lab_id))
    rownames(sp_drug_sen_df) = sp_drug_sen_df$lab_id
    sp_sample = intersect(sp_drug_sen_df$lab_id,colnames(exp_mat))

    #sen_sample = (sp_drug_sen_df[sp_drug_sen_df$lab_id %in% sp_sample,] %>% filter(auc<182.3626))$lab_id ### should change
    #unsen_sample = (sp_drug_sen_df[sp_drug_sen_df$lab_id %in% sp_sample,] %>% filter(auc>232.7781))$lab_id ### should change
    sen_cutoff = sp_drug_sen_df[top_bottom_number,]$auc
    unsen_cutoff = sp_drug_sen_df[
        dim(sp_drug_sen_df)[1]-top_bottom_number,]$auc


    sen_sample = as.vector((sp_drug_sen_df[sp_drug_sen_df$lab_id %in% sp_sample,] %>% filter(auc<sen_cutoff))$lab_id )
    unsen_sample = as.vector((sp_drug_sen_df[sp_drug_sen_df$lab_id %in% sp_sample,] %>% filter(auc>unsen_cutoff))$lab_id )



    df_fpr_plot = data.frame(sample=colnames(exp_mat),exp=t(exp_mat)[,'CEBPA'],drug_sen=rep('NA',dim(exp_mat)[2]))

    df_fpr_plot[sen_sample,]$drug_sen = 'sen'

    df_fpr_plot[unsen_sample,]$drug_sen = 'unsen'

    df_fpr_plot = df_fpr_plot%>%filter(drug_sen!='NA')



    df_fpr_plot$drug_sen = as.factor(df_fpr_plot$drug_sen)

    # p = ggplot(df_fpr_plot,aes(x = drug_sen, y = exp,fill=drug_sen)) + geom_boxplot()
    p <- ggboxplot(df_fpr_plot, x = 'drug_sen', y = 'exp',add = "jitter",color='drug_sen', font.label=) + stat_compare_means()
    p = p + theme(text=element_text(size=18,  family="Arial"))
    p = p + ylab('CEBPA(rpkm)') + xlab(drug)
    plot(p)
    return(p)
}


drug = our_drugs[3]

p = jn_boxplot(drug,exp_mat,drug_sen_df)
















