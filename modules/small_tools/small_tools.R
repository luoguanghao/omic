

if(FALSE){ # flt3 or not ~ crispr ceres
    crispr_df = read_tsv('/home/lgh/my_project/data/depmap/cleaned/CRISPR_gene_effect_cleaned.tsv')
    mut_df = read_tsv('/home/lgh/my_project/data/depmap/cleaned/CCLE_mutations_cleaned.tsv')

    df_for_plot[df_for_plot$sample%in%(mut_df%>%filter(Hugo_Symbol=='FLT3'))$sample,]$flt3 = 'flt3'
    df_for_plot$CHEK1 = as.numeric(df_for_plot$CHEK1)

    ggviolin(df_for_plot,x='flt3',y='CHEK1',
            add='jitter',color='flt3',
            ylab='chek1_ceres')+stat_compare_means()
}

if(FALSE){ # itd or not ~ expression
    exp_df = read_tsv('/home/lgh/my_project/data/VIZOME/cleaned/VIZOME_log2fpkm_expression_cleaned.tsv')
    mut_df = read_tsv('/home/lgh/my_project/data/VIZOME/cleaned/VIZOME_mutation.maf')

    df_for_plot[df_for_plot$sample%in%(mut_df%>%filter(Hugo_Symbol=='FLT3'))$sample,]$flt3 = 'flt3'
    df_for_plot$CHEK1 = as.numeric(df_for_plot$CHEK1)

    ggviolin(df_for_plot,x='flt3',y='CHEK1',
            add='jitter',color='flt3',
            ylab='chek1_ceres')+stat_compare_means()
}


## 做相关，mut与表达或者crispr score做相关
# exp_or_crispr_df是df,第一列是Hugo_Symbol，后面是样本的值：用CRISPR_gene_effect_cleaned_2.tsv
# gene 给exp_or_crispr_df的基因名
# mut_gene 给突变的基因名，也可以给itd
mut_exp_or_crispr <- function(exp_or_crispr_df,mut_df,gene,mut_gene='itd'){
    df_for_sel_sample = data.frame(t(exp_or_crispr_df%>%column_to_rownames('Hugo_Symbol')))%>%rownames_to_column('sample')

    df_for_plot = data.frame(sample=unique(mut_df$sample), mut='WT')

    df_for_plot = merge(df_for_plot,df_for_sel_sample[,c('sample',gene)])
    colnames(df_for_plot)[3] = 'crispr_or_exp'

    if(mut_gene=='itd'){
        df_for_plot[df_for_plot$sample%in%(mut_df%>%filter(ITDorNOT=='ITD'))$sample,]$mut = 'FLT3ITD'
    }else{
        df_for_plot[df_for_plot$sample%in%(mut_df%>%filter(Hugo_Symbol==mut_gene))$sample,]$mut = 'mut'
    }
    
    df_for_plot$crispr_or_exp = as.numeric(df_for_plot$crispr_or_exp)

    ggviolin(df_for_plot, x='mut', y='crispr_or_exp',
            add='jitter', color='mut',
            xlab=mut_gene, ylab=gene)+stat_compare_means()  
}




# 获得一个df，给出哪些样本时突变哪些是野生
# 获取共突变
## 根据这个函数，把前面的两个函数改了！！！！！！！！！！
get_mutation_wt_sample <- function(mut_df,gene){
    if(length(gene)==2){
        mut_stat_df = data.frame(sample=unique(mut_df$Tumor_Sample_Barcode),mutORwt=sprintf('%s_WT',paste(gene,collapse='_')))
        mut_stat_df[mut_stat_df$sample%in%
                unique((mut_df%>%filter(Hugo_Symbol==gene[1]))$Tumor_Sample_Barcode),]$mutORwt=sprintf('%s_mut',gene[1])
        mut_stat_df[mut_stat_df$sample%in%
                unique((mut_df%>%filter(Hugo_Symbol==gene[2]))$Tumor_Sample_Barcode),]$mutORwt=sprintf('%s_mut',gene[2])
        mut_stat_df[mut_stat_df$sample%in%unique((mut_df%>%filter(Hugo_Symbol==gene[1]))$Tumor_Sample_Barcode)&
                    mut_stat_df$sample%in%unique((mut_df%>%filter(Hugo_Symbol==gene[2]))$Tumor_Sample_Barcode),
                ]$mutORwt=sprintf('%s_mut',paste(gene,collapse='_'))
    }else{
        mut_stat_df = data.frame(sample=unique(mut_df$Tumor_Sample_Barcode),mutORwt=sprintf('%s_WT',gene))
        # return(mut_stat_df)
        mut_stat_df[mut_stat_df$sample%in%
                unique((mut_df%>%filter(Hugo_Symbol==gene))$Tumor_Sample_Barcode),]$mutORwt=sprintf('%s_mut',gene)
    }
    return(mut_stat_df)
}















