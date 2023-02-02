#
# cor with drug sen
# cor with survival
source('/home/lgh/my_project/omic/modules/plot/survival.R')
library(tidyverse)
library(broom)
#####################################
#####################################
# cor with survival
# 一个个基因与生存做相关  生存分析
# !!!!!!!!!!!! exp_mat 这里要改
#
# suv_df
# exp_mat： gene*sample
# gene_vec
#
survival_gene_cor <- function(suv_df,exp_df,gene_vec=NULL){

    if(is.null(gene_vec)){
        gene_vec = exp_df$Hugo_Symbol
    }

    print(length(gene_vec))

    df_for_suv_cal = merge(suv_df[,c('sample','OS','OSS')],
        data.frame(t(exp_df%>%filter(Hugo_Symbol%in%gene_vec)%>%
            column_to_rownames('Hugo_Symbol')))%>%rownames_to_column('sample'))


    res_df = data.frame(
        genes=colnames(df_for_suv_cal)[4:length(colnames(df_for_suv_cal))],
        high_exp=88,
        low_exp=88,
        high_low_fc=88,
        p=88
        )

    for(i in 4:length(colnames(df_for_suv_cal))){
        
        sur_res = try(  plot_survival_contin(df_for_suv_cal[,
                c('sample','OS','OSS',colnames(df_for_suv_cal)[i])], colnames(df_for_suv_cal)[i])  )
        
        if("try-error" %in% class(sur_res)){
            next
        }else{
            res_df[i-3, 'p']=sur_res$pval
            res_df[i-3, 'high_exp']=sur_res$surv_diff$exp[1]
            res_df[i-3, 'low_exp']=sur_res$surv_diff$exp[2]
            res_df[i-3, 'high_low_fc']=sur_res$surv_diff$exp[1]/sur_res$surv_diff$exp[2]
        }
    }

    res_df = res_df[order(res_df$p),]
    raw_res_df = res_df
    res_df = res_df%>%filter(p!=88)
    res_df$p.adj = stats::p.adjust(res_df$p, method = "fdr")

    return(list(res_df=res_df,raw_res_df=raw_res_df,df_for_suv_cal=df_for_suv_cal))
}

survival_gene_plot <- function(survival_gene_res_df,df_for_suv_cal,n){
    for(i in 1:n){
        g = survival_gene_res_df$genes[i]
        res = plot_survival_contin(df_for_suv_cal[,c('sample','OS','OSS',g)], g, plot=TRUE)
    }

}




survival_gene_coxph <- function(suv_df,exp_df,gene_vec=NULL,gene_sample=TRUE,OS_RFS='OS'){
    # if(!is.null(gene_vec)) suv_df = suv_df%>%filter(sample%in%)
    if(OS_RFS=='OS'){
        if(gene_sample==TRUE){
            df_for_suv_cal = merge(suv_df[,c('sample','OS','OSS')],
                data.frame(t(exp_df%>%
                    column_to_rownames('Hugo_Symbol')))%>%rownames_to_column('sample'))
        }else{
            df_for_suv_cal = merge(suv_df[,c('sample','OS','OSS')], exp_df)        
        }

        df_for_suv_cal = df_for_suv_cal%>%filter()

        covariates <- colnames(df_for_suv_cal)[4:ncol(df_for_suv_cal)]

        univ_formulas <- sapply(covariates,
                                function(x) as.formula(paste('Surv(OS, OSS)~', x)))
    }else{
        if(gene_sample==TRUE){
            df_for_suv_cal = merge(suv_df[,c('sample','RFS','RFSS')],
                data.frame(t(exp_df%>%
                    column_to_rownames('Hugo_Symbol')))%>%rownames_to_column('sample'))
        }else{
            df_for_suv_cal = merge(suv_df[,c('sample','RFS','RFSS')], exp_df)        
        }

        df_for_suv_cal = df_for_suv_cal%>%filter()

        covariates <- colnames(df_for_suv_cal)[4:ncol(df_for_suv_cal)]

        univ_formulas <- sapply(covariates,
                                function(x) as.formula(paste('Surv(RFS, RFSS)~', x)))
    }

    univ_models <- lapply( univ_formulas, 
        function(x){
            coxph_res = try( coxph(x, data = df_for_suv_cal) )

            if("try-error" %in% class(coxph_res)){
                '-'
            }else{        
                coxph_res
            }
        }
    )


    univ_results <- lapply(univ_models,
                        function(x){
                            if(x=='-'){
                                res<-c(88, 88, 88, 88)
                                names(res)<-c("beta", "HR (95% CI for HR)","wald.test", "p.value")
                            }else{
                                x <- summary(x)
                                p.value<-signif(x$wald["pvalue"], digits=2)
                                wald.test<-signif(x$wald["test"], digits=2)
                                beta<-signif(x$coef[1], digits=2);#coeficient beta
                                HR <-signif(x$coef[2], digits=2);#exp(beta)
                                HR.confint.lower <- signif(x$conf.int[,"lower .95"],2)
                                HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                                HR <- paste0(HR, " (", 
                                            HR.confint.lower, "-", HR.confint.upper, ")")
                                res<-c(beta, HR, wald.test, p.value)
                                names(res)<-c("beta", "HR (95% CI for HR)","wald.test", "p.value")
                            }
                            
                            return(res)
                            #return(exp(cbind(coef(x),confint(x))))
                            })

    res <- t(as.data.frame(univ_results, check.names = FALSE))
    res = as.data.frame(res)%>%rownames_to_column('gene')

    raw_res = res

    res = res[order(res$p.value),]
    res = res%>%filter(p.value!=88)
    res = na.omit(res)
    res$p.adj = stats::p.adjust(res$p.value, method = "fdr")

    count_df = data.frame(count=apply(df_for_suv_cal[,c(4:ncol(df_for_suv_cal))],2,function(x) sum(!is.na(x))))%>%rownames_to_column('gene')
    res = merge(res,count_df,by='gene')

    res = res[order(res$p.value),]
    return(list(res=res, raw_res=raw_res))

}





#####################################
#####################################
# cor with drug sen
# if is gene*sample set t=TRUE
# 一个个基因与药敏做相关
drugsen_gene_cor <- function(drug_sen_df,exp_df,drug,genes=NULL,method=method,aucic50='ic50',t=TRUE,orderby='p'){
    if(t){
        exp_df_t = data.frame(t(exp_df[,-1]))
        colnames(exp_df_t) = exp_df[[1]]
        exp_df_t = exp_df_t%>%rownames_to_column('sample')
        exp_df = exp_df_t
    }
    if(!is.null(genes)){
        exp_df = exp_df[,c('sample',intersect(colnames(exp_df)[2:ncol(exp_df)],genes))]
    }

    df_for_cal = merge((drug_sen_df%>%filter(inhibitor==drug))[,c('sample','inhibitor',aucic50)],exp_df,by='sample')

    drugsen_gene_corr_df = list(drug=c(),gene=c(),r=c(),p=c())
    for(i in 4:ncol(df_for_cal)){
        if(dim(na.omit(df_for_cal[,c(3,i)]))[1]<3){
            next;
        }
        #print(df_for_cal[[i]])
        #print(df_for_cal[[3]])

        cor_res = cor.test(df_for_cal[[3]],df_for_cal[[i]],method=method)
        drugsen_gene_corr_df$drug = c(drugsen_gene_corr_df$drug, drug)
        drugsen_gene_corr_df$gene = c(drugsen_gene_corr_df$gene, colnames(df_for_cal)[i])
        drugsen_gene_corr_df$r = c(drugsen_gene_corr_df$r, cor_res$estimate[[1]])
        drugsen_gene_corr_df$p = c(drugsen_gene_corr_df$p, cor_res$p.value)
    }
    drugsen_gene_corr_df = data.frame(drugsen_gene_corr_df)
    drugsen_gene_corr_df$p.adj = p.adjust(drugsen_gene_corr_df$p, method = 'fdr')
    drugsen_gene_corr_df = drugsen_gene_corr_df[order(drugsen_gene_corr_df$p),]

    # ===== add count

    # df_for_cal = all_res_ls[['367C']]$phosph$df_for_cal
    rownames(df_for_cal) = NULL
    drugsen_gene_corr_df_with_count = drugsen_gene_corr_df
    drugsen_gene_corr_df_with_count = merge(drugsen_gene_corr_df_with_count,
                            data.frame(count=apply(df_for_cal[,c(4:ncol(df_for_cal))],2,function(x) sum(!is.na(x))))%>%rownames_to_column('gene'),
                            by='gene')

    # ===== add group
    drugsen_gene_corr_df_with_count_with_group = rbind(drugsen_gene_corr_df_with_count['NU',],
                                drugsen_gene_corr_df_with_count)
    drugsen_gene_corr_df_with_count_with_group[1,] = aucic50

    drugsen_gene_corr_df_with_count_with_group$gene = factor(drugsen_gene_corr_df_with_count_with_group$gene,
                                    level=drugsen_gene_corr_df_with_count_with_group$gene)
    # ===== add sample matrix
    drugsen_gene_corr_df_with_count_with_group = merge(
        drugsen_gene_corr_df_with_count_with_group,
        ( data.frame(t(df_for_cal%>%column_to_rownames('sample')))%>%rownames_to_column('gene') )[-1,],
        by='gene'
    )
    drugsen_gene_corr_df_with_count_with_group = 
        drugsen_gene_corr_df_with_count_with_group[order(drugsen_gene_corr_df_with_count_with_group$gene),]
    # test_df_with_count = test_df_with_count

    if(orderby=='p'){
        tmp_df = drugsen_gene_corr_df_with_count_with_group
        drugsen_gene_corr_df_with_count_with_group = rbind(tmp_df[1,],tmp_df[2:nrow(tmp_df),][order(tmp_df[2:nrow(tmp_df),]$p),])
    }
    # =====

    return(list(drugsen_gene_corr_df=drugsen_gene_corr_df_with_count,
                drugsen_gene_corr_df_with_count=drugsen_gene_corr_df_with_count_with_group,
                df_for_cal=df_for_cal))
}


# cor with drug sen
# 一个个基因与连续变量表型做相关
XXX.drugsen_gene_cor <- function(drug_sen_df,exp_df,drug,genes=NULL,method='spearman',aucic50='ic50',t=TRUE){
    if(t){
        exp_df_t = data.frame(t(exp_df[,-1]))
        colnames(exp_df_t) = exp_df[[1]]
        exp_df_t = exp_df_t%>%rownames_to_column('sample')
        exp_df = exp_df_t
    }
    if(!is.null(genes)){
        exp_df = exp_df[,c('sample',intersect(colnames(exp_df)[2:ncol(exp_df)],genes))]
    }

    df_for_cal = merge((drug_sen_df%>%filter(inhibitor==drug))[,c('sample','inhibitor','ic50')],exp_df,by='sample')    

    drugsen_gene_corr_df = list(drug=c(),gene=c(),r=c(),p=c())
    for(i in 4:ncol(df_for_cal)){
        if(length(na.omit(df_for_cal[[i]]))<3){
            next;
        }
        #print(df_for_cal[[i]])
        #print(df_for_cal[[3]])

        cor_res = cor.test(df_for_cal[[3]],df_for_cal[[i]],method='spearman')
        drugsen_gene_corr_df$drug = c(drugsen_gene_corr_df$drug, drug)
        drugsen_gene_corr_df$gene = c(drugsen_gene_corr_df$gene, colnames(df_for_cal)[i])
        drugsen_gene_corr_df$r = c(drugsen_gene_corr_df$r, cor_res$estimate[[1]])
        drugsen_gene_corr_df$p = c(drugsen_gene_corr_df$p, cor_res$p.value)
    }
    drugsen_gene_corr_df = data.frame(drugsen_gene_corr_df)
    drugsen_gene_corr_df$p.adj = p.adjust(drugsen_gene_corr_df$p, method = 'fdr')
    drugsen_gene_corr_df = drugsen_gene_corr_df[order(drugsen_gene_corr_df$p),]
    return(list(drugsen_gene_corr_df=drugsen_gene_corr_df,
                df_for_cal=df_for_cal))
}




#####################################
#####################################
# cor with drug sen
# 一个个基因突变与生存做相关
survival_mut_test <- function(){

}


# 一个个基因突变与药敏做 ttest
# mut_mat_df : c('sample','gene1','gene2')
drugsen_mut_test <- function(drug_sen_df,mut_mat_df,drug,genes=NULL,method='spearman',aucic50='ic50'){
    
    df_for_cal = merge((drug_sen_df%>%filter(inhibitor==drug))[,c('sample','inhibitor','ic50')],mut_mat_df,by='sample')
    
    drugsen_mut_test_df = list(drug=c(),gene=c(),fc_median=c(),p=c())
    for(i in 4:ncol(df_for_cal)){
        if(length(unique(df_for_cal[[i]]))!=1){
            colnames(df_for_cal)[i]
            test_res = wilcox.test(df_for_cal[df_for_cal[[i]]==1,]$ic50,
                                   df_for_cal[df_for_cal[[i]]==0,]$ic50)
            median_mut = median(df_for_cal[df_for_cal[[i]]==1,]$ic50)
            median_wt = median(df_for_cal[df_for_cal[[i]]==0,]$ic50)
            drugsen_mut_test_df$drug = c(drugsen_mut_test_df$drug, drug)
            drugsen_mut_test_df$gene = c(drugsen_mut_test_df$gene, colnames(df_for_cal)[i])
            drugsen_mut_test_df$fc_median = c(drugsen_mut_test_df$fc_median, log(median_mut/median_wt))
            drugsen_mut_test_df$p = c(drugsen_mut_test_df$p, test_res$p.value)
        }
    }
    drugsen_mut_test_df = data.frame(drugsen_mut_test_df)
    drugsen_mut_test_df$p.adj = p.adjust(drugsen_mut_test_df$p, method = 'fdr')
    drugsen_mut_test_df = drugsen_mut_test_df[order(drugsen_mut_test_df$p),]
    return(list(drugsen_mut_test_df=drugsen_mut_test_df,
                df_for_cal=df_for_cal))    
}




#####################################
#####################################
# diff exp onr by one
# exp_df: sample*gene   ###
# group_df: sample,group
# t: tran gene*sample into sample*gene
# method: wilcox, 后期加上t.text
#
# 加入compare_means的快速方法
# ===
# df_for_cal: sample,group,gene1,gene2...
# 比较: g2 vs g1
diff_exp_one_by_one <- function(exp_df, group_df, g1,g2, 
                method='wilcox', ggpubr=TRUE, genes=NULL, 
                t=FALSE, log=FALSE){

    if(t){
        exp_df_t = data.frame(t(exp_df[,-1]))
        colnames(exp_df_t) = exp_df[[1]]
        exp_df_t = exp_df_t%>%rownames_to_column('sample')
        exp_df = exp_df_t
    }
    if(!is.null(genes)){
        exp_df = exp_df[,c('sample',intersect(colnames(exp_df)[2:ncol(exp_df)],genes))]
    }

    df_for_cal = merge(group_df,exp_df,by='sample') # 只剩有用的!!!

    test_df = list(gene=c(),fc_median=c(),p=c())

    print(sprintf('compare: g2/g1: %s/%s',g2,g1))

    if(!ggpubr){
        for(i in 3:ncol(df_for_cal)){

            if(length(na.omit(df_for_cal[[i]]))<3){
                next;
            }
            
            test_res = wilcox.test(df_for_cal[df_for_cal$group==g1,][[i]],
                                    df_for_cal[df_for_cal$group==g2,][[i]],exact=TRUE)
            median_g1 = median(df_for_cal[df_for_cal$group==g1,][[i]])
            median_g2 = median(df_for_cal[df_for_cal$group==g2,][[i]])

            test_df$gene = c(test_df$gene, colnames(df_for_cal)[i])

            if(log==TRUE){
                test_df$fc_median = c(test_df$fc_median, 2**(median_g2-median_g1))
            }else{
                test_df$fc_median = c(test_df$fc_median, median_g2/median_g1)
            }
            test_df$log2_fc = log2(test_df$fc_median)
            test_df$p = c(test_df$p, test_res$p.value)
        }

        test_df = data.frame(test_df)
        test_df$p.adj = p.adjust(test_df$p, method = 'fdr')
        test_df = test_df[order(test_df$p),]
        test_df = test_df[order(test_df$p),c('gene','log2_fc','p','p.adj')]

    }else{
        
        data3 <- reshape2::melt(df_for_cal, 
                        id.vars = c("sample","group"))
        colnames(data3) = c('sample','group','gene','value')
        data3 = na.omit(data3)
        tt=data3%>%group_by(group,gene)%>%summarise(n=n())
        paichu = as.vector((tt%>%filter(gene%in%names(table(tt$gene))[table(tt$gene)!=2]))$gene)
        data3 = data3%>%filter(!(gene%in%paichu))

        # 要排除掉不行的。。。

        test_res = ggpubr::compare_means(value~group,data3,group.by='gene',exact=TRUE)

        g1_median = data3%>%filter(group==g1) %>%
        group_by(gene) %>%
        summarise(median = median(value))

        g2_median = data3%>%filter(group==g2) %>%
        group_by(gene) %>%
        summarise(median = median(value))
        #View(head(data3))
        #View(head(g1_median))
        #View(head(g2_median))

        
        #View(head(test_df))

        if(log==TRUE){
            test_df = merge(test_res[,c('gene','p','p.adj')],
                            data.frame(gene=g1_median$gene,
                                    fc_median = 2**(g2_median$median-g1_median$median)))
        }else{
            test_df = merge(test_res[,c('gene','p','p.adj')],
                            data.frame(gene=g1_median$gene,
                                    fc_median = g2_median$median/g1_median$median))

        }
        test_df = test_df[order(test_df$p),c('gene','fc_median','p','p.adj')]
        test_df$log2_fc = log2(test_df$fc_median)
        test_df = test_df[order(test_df$p),c('gene','log2_fc','p','p.adj')]

    }


    # === add count

    data4 <- reshape2::melt(df_for_cal, 
        id.vars = c("sample",'group'))
    data4 = na.omit(data4)

    data4_count = data4%>%
        group_by(group,variable) %>%
        summarise(count = n())

    data5 = reshape2::dcast(data4_count, variable~group)[,c('variable',g1,g2)]

    #data5[is.na(data5[['RUNX1:RUNX1T1']]),][['RUNX1:RUNX1T1']]=0

    test_df_with_count = merge(test_df,
        data5,by.x='gene',by.y='variable')
    test_df_with_count$gene = as.character(
        test_df_with_count$gene
    )
    
    test_df_with_count = test_df_with_count[order(test_df_with_count$log2_fc,decreasing=TRUE),]

    # === add group
    
    test_df_with_count_with_group = rbind(test_df_with_count['NU',],
                                test_df_with_count)
    test_df_with_count_with_group[1,] = 'group'

    test_df_with_count_with_group$gene = factor(test_df_with_count_with_group$gene,
                                    level=test_df_with_count_with_group$gene)

    # --- add mat

    df_for_cal = df_for_cal[order(df_for_cal$group),]
    rownames(df_for_cal) = NULL

    tmp_df = merge(
        test_df_with_count_with_group,
        data.frame(t(df_for_cal%>%column_to_rownames('sample')))%>%rownames_to_column('gene'),
        by='gene'
    )
    tmp_df = tmp_df[order(tmp_df$gene),]
    test_df_with_count_with_group = tmp_df

    # ===


    return(list(test_df=test_df_with_count, # stat value; group count; matrix
                test_df_with_count=test_df_with_count_with_group, # stat value; group count; matrix; row of group
                df_for_cal=df_for_cal))  

}


#####################################
#####################################
# kw exp onr by one
# exp_df: sample*gene   ###
# group_df: sample,group
# t: tran gene*sample into sample*gene
## if df is gene*sample, set t=TRUE
# method: wilcox, 后期加上t.text
#
# 加入compare_means的快速方法
# ===
# df_for_cal: sample,group,gene1,gene2...
#

kw_exp_one_by_one <- function(exp_df, group_df, 
                method='wilcox', ggpubr=TRUE, genes=NULL, 
                t=FALSE, log=FALSE,median_mean='median'){
    # protein_df
    if(t){
        #exp_df_t = data.frame(t(exp_df[,-1]))
        #colnames(exp_df_t) = exp_df[[1]]
        #exp_df_t = exp_df_t%>%rownames_to_column('sample')
        #exp_df = exp_df_t
        colnames(exp_df)[1] = 'gene'
        data3 <- reshape2::melt(exp_df, id.vars = c("gene"))
        colnames(data3) = c('variable','sample','value')
    }else{
        data3 <- reshape2::melt(exp_df, id.vars = c("sample"))

    }
    if(!is.null(genes)){
        # exp_df = exp_df[,c('sample',intersect(colnames(exp_df)[2:ncol(exp_df)],genes))]
    }

    # data3 <- reshape2::melt(exp_df, id.vars = c("sample"))

    df_for_cal_pro = merge(group_df,data3)%>%filter(!is.na(value))
    # with pvalues and further stats
    kw_res_df_pro = df_for_cal_pro%>%filter() %>% 
    nest(-variable) %>% 
    mutate(kw=map(data,~kruskal.test(.x$value ~ .x$group))) %>%
    mutate(tidied = map(kw, tidy)) %>% 
    unnest(tidied, .drop = T)

    kw_res_df_pro$p.adj = p.adjust(kw_res_df_pro$p.value)
    #group_median_pro = df_for_cal_pro%>%group_by(variable,group)%>%summarize(median=median(value))
    #write_tsv(group_median_pro,file=sprintf('clustering_analysis_data_res/xjy_kw_analysis/%s_group_median_pro.tsv',param_ls[i]))
    #write_tsv(merge(kw_res_df_pro[,c(-2,-3)],reshape2::dcast(group_median_pro, variable~group, value.var='median'),by='variable'),
    #          file=sprintf('%s/clustering_analysis_data_res/xjy_kw_analysis/%s_kw_res_df_pro.tsv',cluster_res_dir,param_ls[i]))
    group_count_pro = df_for_cal_pro%>%group_by(variable,group)%>%summarize(count=n())
    tmp_count = reshape2::dcast(group_count_pro, variable~group, value.var='count')
    colnames(tmp_count)[2:4] = paste('count_',colnames(tmp_count)[2:4],sep='')

    if(median_mean=='median'){
        group_median_pro = df_for_cal_pro%>%group_by(variable,group)%>%summarize(median=median(value))
        tmp_median = reshape2::dcast(group_median_pro, variable~group, value.var='median')
        colnames(tmp_median)[2:4] = paste('median_',colnames(tmp_median)[2:4],sep='')

        tmp_1 = merge(kw_res_df_pro[,c(-2,-3)],tmp_median,by='variable')
        tmp_1 = merge(tmp_1,tmp_count,by='variable')
        colnames(tmp_1)[1] = 'gene'    
    }else{
        group_mean_pro = df_for_cal_pro%>%group_by(variable,group)%>%summarize(mean=mean(value))
        tmp_mean = reshape2::dcast(group_mean_pro, variable~group, value.var='mean')
        colnames(tmp_mean)[2:4] = paste('mean_',colnames(tmp_mean)[2:4],sep='')

        tmp_1 = merge(kw_res_df_pro[,c(-2,-3)],tmp_mean,by='variable')
        tmp_1 = merge(tmp_1,tmp_count,by='variable')
        colnames(tmp_1)[1] = 'gene' 
    }

    if(t){
        tmp_mat = data.frame(exp_df)[,c('gene',group_df$sample)]
    }else{
        tmp_mat = data.frame(t(exp_df%>%column_to_rownames('sample')))[,group_df$sample]
        # colnames(tmp_mat) = paste(group_df$group,colnames(tmp_mat),sep='_')
        tmp_mat = tmp_mat%>%rownames_to_column('gene')      
    }
   
    tmp_1_with_mat = merge(tmp_1,tmp_mat,by='gene')  
    tmp_1_with_mat$gene = as.character(tmp_1_with_mat$gene)
    tmp_1$gene = as.character(tmp_1$gene)
    # tmp_2_with_mat = merge(tmp_2,tmp_mat,by='gene')

    # add anno #
    tmp_tmp_1_with_mat = rbind(tmp_1_with_mat['NU',],tmp_1_with_mat)
    tmp_tmp_1_with_mat[1,] = 'group'
    tmp_tmp_1_with_mat[1,(ncol(tmp_1)+1):ncol(tmp_tmp_1_with_mat)] = 
        (group_df%>%column_to_rownames('sample'))[colnames(tmp_tmp_1_with_mat)[(ncol(tmp_1)+1):ncol(tmp_tmp_1_with_mat)],]
    ##

    #write_tsv(tmp_1_with_mat,
    #        file=sprintf('ELN_analysis_result/eln_kw/kw_res_df_pro.tsv'))
    #write_tsv(tmp_tmp_1_with_mat,
    #        file=sprintf('ELN_analysis_result/eln_kw/kw_res_df_pro_anno_group.tsv'))


    tmp_1 = tmp_1[order(tmp_1$p.value),]
    tmp_1$p.adj = p.adjust(tmp_1$p.value)
    return(list(
        kw_df=tmp_1, # stat value; group count; matrix
        kw_df_anno_group=tmp_tmp_1_with_mat, # stat value; group count; matrix; row of group
        df_for_cal=df_for_cal_pro
    ))

}

if(FALSE){
    kinase_anno = data.frame(gene=paste(KSEAapp::KSData[,c('SUB_GENE')],KSEAapp::KSData[,c('SUB_MOD_RSD')],sep='_'),
                            kinase=as.character(KSEAapp::KSData$KINASE),
                            kinase_gene=as.character(KSEAapp::KSData$GENE),
                            kinase_id=as.character(KSEAapp::KSData$KIN_ACC_ID)
                        )
    group_df = eln_df%>%filter(ELN!='ELNNA') # !!!!!!
    colnames(group_df)[2] = 'group'
}
kw_phosphsite_one_by_one <- function(phosph_df, group_df, kinase_anno, 
                method='wilcox', ggpubr=TRUE, genes=NULL, 
                t=FALSE, log=FALSE){

    if(t){
        exp_df_t = data.frame(t(phosph_df[,-1]))
        colnames(exp_df_t) = exp_df[[1]]
        exp_df_t = exp_df_t%>%rownames_to_column('sample')
        phosph_df = exp_df_t
    }
    if(!is.null(genes)){
        exp_df = exp_df[,c('sample',intersect(colnames(exp_df)[2:ncol(exp_df)],genes))]
    }

    data3 <- reshape2::melt(phosph_df, id.vars = c("sample"))
    df_for_cal_phosph_site = merge(group_df,data3)%>%filter(!is.na(value))

    kw_res_df_phosph_site = df_for_cal_phosph_site%>%filter() %>% 
    nest(-variable) %>% 
    mutate(kw=map(data,~kruskal.test(.x$value ~ .x$group))) %>%
    mutate(tidied = map(kw, tidy)) %>% 
    unnest(tidied, .drop = T)

    ############
    # xu output
    group_count_phosph_site = df_for_cal_phosph_site%>%group_by(variable,group)%>%summarize(count=n())
    tmp_count = reshape2::dcast(group_count_phosph_site, variable~group, value.var='count')
    colnames(tmp_count)[2:4] = paste('count_',colnames(tmp_count)[2:4])

    group_median_phosph_site = df_for_cal_phosph_site%>%group_by(variable,group)%>%summarize(median=median(value))
    tmp_median = reshape2::dcast(group_median_phosph_site, variable~group, value.var='median')
    colnames(tmp_median)[2:4] = paste('median_',colnames(tmp_median)[2:4])

    tmp_1 = merge(kw_res_df_phosph_site[,c(-2,-3)],tmp_median,by='variable')
    tmp_1 = merge(tmp_1,tmp_count,by='variable')
    colnames(tmp_1)[1] = 'gene'    
    tmp_2 = merge(kinase_anno,tmp_1,by.x='gene',by.y='gene',all.y=TRUE)    

    tmp_mat = data.frame(t(phosph_df%>%column_to_rownames('sample')))[,group_df$sample]
    # colnames(tmp_mat) = paste(group_df$group,colnames(tmp_mat),sep='_')
    tmp_mat = tmp_mat%>%rownames_to_column('gene')    

    tmp_1_with_mat = merge(tmp_1,tmp_mat,by='gene')    
    tmp_2_with_mat = merge(tmp_2,tmp_mat,by='gene')
    tmp_1_with_mat$gene = as.character(tmp_1_with_mat$gene)
    tmp_2_with_mat$gene = as.character(tmp_2_with_mat$gene)
    tmp_1$gene = as.character(tmp_1$gene)
    tmp_2$gene = as.character(tmp_2$gene)

    # add sample anno #

    tmp_tmp_1_with_mat = rbind(tmp_1_with_mat['NU',],tmp_1_with_mat)
    tmp_tmp_1_with_mat[1,] = 'group'
    tmp_tmp_1_with_mat[1,(ncol(tmp_1)+1):ncol(tmp_tmp_1_with_mat)] = 
        (group_df%>%column_to_rownames('sample'))[colnames(tmp_tmp_1_with_mat)[(ncol(tmp_1)+1):ncol(tmp_tmp_1_with_mat)],]

    tmp_tmp_2_with_mat = rbind(tmp_2_with_mat['NU',],tmp_2_with_mat)
    tmp_tmp_2_with_mat[1,] = 'group'
    tmp_tmp_2_with_mat[1,(ncol(tmp_2)+1):ncol(tmp_tmp_2_with_mat)] = 
        (group_df%>%column_to_rownames('sample'))[colnames(tmp_tmp_2_with_mat)[(ncol(tmp_2)+1):ncol(tmp_tmp_2_with_mat)],]
    ####
    #write_tsv(tmp_1_with_mat,
    #        file=sprintf('ELN_analysis_result/eln_kw/kw_res_df_phosph_site.tsv'))
    #write_tsv(tmp_2_with_mat,
    #        file=sprintf('ELN_analysis_result/eln_kw/kw_res_df_phosph_site_anno_kinase.tsv'))

    #write_tsv(tmp_tmp_1_with_mat,
    #        file=sprintf('ELN_analysis_result/eln_kw/kw_res_df_phosph_site_anno_group.tsv'))
    #write_tsv(tmp_tmp_2_with_mat,
    #        file=sprintf('ELN_analysis_result/eln_kw/kw_res_df_phosph_site_anno_kinase_anno_group.tsv'))


    return(list(
        kw_df=tmp_1, # stat value; group count; matrix
        kw_df_anno_group=tmp_tmp_1_with_mat, # stat value; group count; matrix; row of group
        kw_df_anno_kinase=tmp_2, # stat value; group count; matrix; anno kinase
        kw_df_anno_kinase_anno_group=tmp_tmp_2_with_mat, # stat value; group count; matrix; row of group; anno kinase
        df_for_cal=df_for_cal_phosph_site
    ))

}


#####################################
#####################################
# fisher one by one










# demo_pipeline
if(FALSE){
    
}



#####################################
#####################################
# anno kinase


anno_kinase <- function(test_df,test_df_with_count){

    kinase_anno = data.frame(gene=paste(KSEAapp::KSData[,c('SUB_GENE')],
                                KSEAapp::KSData[,c('SUB_MOD_RSD')],sep='_'),
                            kinase=KSEAapp::KSData$KINASE,
                            kinase_gene=KSEAapp::KSData$GENE,
                            kinase_id=KSEAapp::KSData$KIN_ACC_ID)
    # ===
    test_df_anno_kinase = merge(kinase_anno,test_df,all.y=TRUE)

    # ===
    kinase_anno_addgroup = rbind(data.frame(gene='group',kinase='group',kinase_gene='group',kinase_id='group'),
        kinase_anno)

    test_df_with_count_anno_kinase = merge(kinase_anno_addgroup,
        test_df_with_count,
        all.y=TRUE)
    test_df_with_count_anno_kinase = rbind(test_df_with_count_anno_kinase%>%filter(gene=='group'),
        test_df_with_count_anno_kinase%>%filter(gene!='group'))
    #write_tsv(tmp_df,
    #          file='survival_1year_2year_res/survival_1year_2year_phosph_test_df_anno_kinase_with_count.tsv')

    return(list(
        test_df_anno_kinase = test_df_anno_kinase,
        test_df_with_count_anno_kinase = test_df_with_count_anno_kinase
    ))
}
