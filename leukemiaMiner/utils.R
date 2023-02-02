####################################
####################################
# Expression Analysis & GSEA
####################################
####################################
# omic_df
# score_df:first column is sample names

get_score_gene_corr_df = function(score_df, omic_df, method='spearman'){

    # df_for_cal = merge((drug_sen_df%>%filter(inhibitor==drug))[,c('sample','inhibitor',aucic50)],exp_df,by='sample')
    df_for_cal = merge(score_df, omic_df, by='sample')

    drugsen_gene_corr_df = list(gene=c(),r=c(),p=c(),count=c())
    for(i in 3:ncol(df_for_cal)){
        if(dim(na.omit(df_for_cal[,c(2,i)]))[1]<3){
            next;
        }
        #print(df_for_cal[[i]])
        #print(df_for_cal[[3]])

        cor_res = cor.test(df_for_cal[[2]],df_for_cal[[i]],method='spearman')
        # drugsen_gene_corr_df$drug = c(drugsen_gene_corr_df$drug, drug)
        drugsen_gene_corr_df$gene = c(drugsen_gene_corr_df$gene, colnames(df_for_cal)[i])
        drugsen_gene_corr_df$r = c(drugsen_gene_corr_df$r, cor_res$estimate[[1]])
        drugsen_gene_corr_df$p = c(drugsen_gene_corr_df$p, cor_res$p.value)
        drugsen_gene_corr_df$count = c(drugsen_gene_corr_df$count, dim(na.omit(df_for_cal[,c(2,i)]))[1])
    }
    drugsen_gene_corr_df = data.frame(drugsen_gene_corr_df)
    drugsen_gene_corr_df$p.adj = p.adjust(drugsen_gene_corr_df$p, method = 'fdr')
    drugsen_gene_corr_df = drugsen_gene_corr_df[order(drugsen_gene_corr_df$p),]
    drugsen_gene_corr_df$method = 'spearman'
    # =====================
    
    

    
    return(drugsen_gene_corr_df)
    
}


### do gsea preranked

if(FALSE){

    PATH_DATA = '/home/lgh/my_project/data'
    gmtfile = list(
        #gobp=paste(PATH_DATA,'/pathway/c5.go.bp.v7.5.1.symbols.gmt',sep=''),
        #kegg=paste(PATH_DATA,'/pathway/c2.cp.kegg.v7.5.1.symbols.gmt',sep=''),
        reactome=paste(PATH_DATA,'/pathway/c2.cp.reactome.v7.5.1.symbols.gmt',sep=''),
        hallmark='/home/lgh/my_project/data/pathway/h.all.v7.5.1.symbols.gmt'#,
        #tft=paste(PATH_DATA,'/pathway/c3.tft.v7.5.1.symbols.add_cebpa.gmt',sep='')
    )

    # ========================

    dir.create(sprintf('%s%s/gsea/','metabolism/','PCA_rna'),recursive=T)

    # ========================

    gsea_list_exp = list()
    #omics = 'exp'
    #for(drug in names(omics_cor_analysis_res_ls)){
    #for(drug in 'ABT-263 (Navitoclax)'){
    tmp_df = get_factor_gene_corr_df_rna%>%filter(count>20)
    # tmp_df[,3:6] = sapply(tmp_df[,3:6], as.numeric)
    tmp_df = tmp_df[,1:2]
    tmp_df = na.omit(tmp_df[order(tmp_df$r),])


    fn = sprintf('%s%s/gsea/%s_%s_%s.rnk','metabolism/','PCA_rna','rna','dim.1','corr')
    write_tsv(tmp_df, 
            file=fn,col_names=FALSE)
    # gsea_list[[drug]] = list()
    for(item in names(gmtfile)){
        gsea_list_exp[[item]] = gsea_cmd_prerank(rnk=fn,
                                                    gmt=gmtfile[[item]],
                                                    rpt_label=sprintf('%s_%s_%s_%s','rna','dim.1','corr',item),
                                                    permute=1000,
                                                    out=sprintf('%s%s/gsea','metabolism/','PCA_rna'))
    }

    gsea_res_list_pro = list()
    # message(drug)
    for( term in names(gsea_list_exp) ){

        gsea_res_list_pro[[term]] = 
            system(gsea_list_exp[[term]], intern = TRUE)

    }



}



####################################
####################################
# Drug Sen Analysis
####################################
####################################
#
#
# cor_drug_factor_df_ch_rna = factor_drug_sen_correlation(drug_sen_df=ch_drug_df,factor_df=factor_analysis_res$factor_df,factor_dim='Dim.1')
# cor_drug_factor_df_ch_rna%>%filter(p<0.05)

factor_drug_sen_correlation <- function(drug_sen_df,factor_df,factor_dim){

    factor_df = factor_df[,c('sample',sprintf('%s',factor_dim))]

    drug_cor_df = merge(factor_df,drug_sen_df)
    drug_cor_df = drug_cor_df%>%filter(inhibitor%in%names(table(drug_cor_df$inhibitor)>=4))
    cor_drug_factor_df = list(drug=c(),r=c(),p=c())

    for(drug_tmp in unique(drug_cor_df$inhibitor)){
        try_res = try({
            cor_res= cor.test(as.formula(sprintf('~%s+ic50',factor_dim)),data=drug_cor_df%>%filter(inhibitor==drug_tmp),method='spearman')
            cor_drug_factor_df$drug = c(cor_drug_factor_df$drug, drug_tmp)
            cor_drug_factor_df$r = c(cor_drug_factor_df$r, cor_res$estimate[[1]])
            cor_drug_factor_df$p = c(cor_drug_factor_df$p, cor_res$p.value)
        })

    }

    cor_drug_factor_df = data.frame(cor_drug_factor_df)
    cor_drug_factor_df$adj.p = p.adjust(cor_drug_factor_df$p)

    return(cor_drug_factor_df)
}



utils.do_coxph <- function(df_for_surv=NULL, value_df=NULL, surv_df=NULL, var, OS_RF){
    
    if(is.null(df_for_surv)) df_for_surv = merge(surv_df,value_df,by='sample')
    # disc_result = plot_survival_contin(df_for_suv_plot, sprintf('Dim.%s',ii), plot=T, title=NULL, OS_RF='OS')\
    coxph_res = coxph(as.formula(sprintf('Surv(%s,%sS)~%s',OS_RF,OS_RF,var)), df_for_surv )
    res_df_coxph_1_tmp =  summary(coxph_res)
    p.value<-signif(res_df_coxph_1_tmp$wald["pvalue"], digits=2)
    wald.test<-signif(res_df_coxph_1_tmp$wald["test"], digits=2)
    beta<-signif(res_df_coxph_1_tmp$coef[1], digits=2);#coeficient beta
    HR <-signif(res_df_coxph_1_tmp$coef[2], digits=2);#exp(beta)
    HR.confint.lower <- signif(res_df_coxph_1_tmp$conf.int[,"lower .95"],2)
    HR.confint.upper <- signif(res_df_coxph_1_tmp$conf.int[,"upper .95"],2)
    HR <- paste0(HR, " (", 
                HR.confint.lower, "-", HR.confint.upper, ")")

    return(list(
        coxph_df=data.frame(var=var,
            wald.test=wald.test,
            beta=beta,HR=HR,
            HR.confint.lower=HR.confint.lower,
            HR.confint.upper=HR.confint.upper,
            p.value=p.value),
        coxph_res=coxph_res
    ))
    
}




utils.do_km_survival <- function(df_for_surv=NULL, value_df=NULL, surv_df=NULL, var, OS_RF){

    if(is.null(df_for_surv)) df_for_surv = merge(surv_df,value_df, by='sample')

    km_res = plot_survival_contin(df_for_surv, var, plot=F, title=NULL, OS_RF='OS')

    return(list(
        km_df=data.frame(var=var,
            exp_high=km_res$surv_diff$exp[1],
            exp_low=km_res$surv_diff$exp[2],
            n_high=km_res$surv_diff$n[[1]],
            n_low=km_res$surv_diff$n[[2]],
            p=km_res$pval),
        km_res=km_res
    ))

}

