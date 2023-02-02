############################################
## Group ANALYSIS
############################################ 

############
## omic
# omic_df: sample*gene
# factor_sel: one factor or 'all'
# type: PCA,MOFA
# m: pearson;spearman

############
## Volcano

discovery.group.volcano <- function(cor_res, show_genes=NULL, thfc=0.25, thp=0.05,log2FC='log2FoldChange', adj='padj'){

    df_vol = cor_res

    colnames(df_vol)[1] = 'gene'

    df_vol$p_show = df_vol[[adj]]
    df_vol$log2FC = df_vol[[log2FC]]

    df_vol$change = 'grey'
    try_res = try({df_vol[(df_vol$log2FC< -thfc)&(df_vol[[adj]]< thp),]$change = 'blue'})
    try_res = try({df_vol[(df_vol$log2FC> thfc)&(df_vol[[adj]]< thp),]$change = 'red'})
    df_vol$change = factor(df_vol$change,level=c('blue','grey','red'))
    print(table(df_vol$change))

    # === 
    p <- ggplot(
        # 数据、映射、颜色
        df_vol, aes(x = log2FC, y = -log10(p_show), colour=change)) +
        geom_point(alpha=0.4, size=2) +
        scale_color_manual(values=names(table(df_vol$change)))+
        # scale_color_manual(values=c('blue','grey','red'))+
        # 辅助线
        geom_vline(xintercept=c(-thfc, thfc),lty=4,col="black",lwd=0.8) +
        geom_hline(yintercept = -log10(thp),lty=4,col="black",lwd=0.8) +
        # 坐标轴
        labs(x=log2FC,
            y=sprintf("-log10(%s)",adj))+
        theme_bw()+
        # 图例
        theme(plot.title = element_text(hjust = 0.5), 
                legend.position="right", 
                legend.title = element_blank())

    if(!is.null(show_genes)){
        p = p + ggrepel::geom_label_repel(data = df_vol%>%filter(gene%in%show_genes), 
                        aes(x = log2FC, y = -log10(p_show), label = gene),
                        size = 5, box.padding = unit(0.5, "lines"),
                        point.padding = unit(0.8, "lines"), 
                        segment.color = "black", 
                        show.legend = FALSE) +
            geom_point(data = df_vol%>%filter(gene%in%show_genes),
                    aes(x = log2FC, y = -log10(p_show), label = gene),
                    alpha=0.4, size=3, color='black') +
            theme (legend.position = "none")
    }

    return(p)

}



############
## drug




############################################
############################################
## FACTOR ANALYSIS
############################################ 
############################################
# http://202.127.26.99:6424/notebooks/my_home/my_project/hl/AML/feature_selection/eln_factor_analysis-up_eln_59genes.ipynb

############
## omic
# omic_df: sample*gene
# factor_sel: one factor or 'all'
# type: PCA,MOFA
# m: pearson;spearman
discovery.factor.omic.cor <- function(omic_df, factor_df, factor_sel, type, m='spearman'){
    if(type=='MOFA'){
        df_for_cal = merge((factor_df%>%filter(factor==factor_sel))[,c('sample','value')], omic_df)
    }else if(type=='PCA'){
        df_for_cal = merge(factor_df[,c('sample',factor_sel)],omic_df)
        colnames(df_for_cal)[2] = 'value'
    }else{
        print('ERROR type')
    }

    drugsen_gene_corr_df = list(gene=c(),r=c(),p=c(),count=c())
    for(i in 3:ncol(df_for_cal)){
        if(dim(na.omit(df_for_cal[,c(2,i)]))[1]<3){
            next;
        }
        #print(df_for_cal[[i]])
        #print(df_for_cal[[3]])

        cor_res = cor.test(df_for_cal[[2]],df_for_cal[[i]],method=m)
        # drugsen_gene_corr_df$drug = c(drugsen_gene_corr_df$drug, drug)
        drugsen_gene_corr_df$gene = c(drugsen_gene_corr_df$gene, colnames(df_for_cal)[i])
        drugsen_gene_corr_df$r = c(drugsen_gene_corr_df$r, cor_res$estimate[[1]])
        drugsen_gene_corr_df$p = c(drugsen_gene_corr_df$p, cor_res$p.value)
        drugsen_gene_corr_df$count = c(drugsen_gene_corr_df$count, dim(na.omit(df_for_cal[,c(2,i)]))[1])
    }
    drugsen_gene_corr_df = data.frame(drugsen_gene_corr_df)
    drugsen_gene_corr_df$p.adj = p.adjust(drugsen_gene_corr_df$p, method = 'fdr')
    drugsen_gene_corr_df = drugsen_gene_corr_df[order(drugsen_gene_corr_df$p),]
    drugsen_gene_corr_df$method = m
    # =====================
    
    drugsen_gene_corr_df = drugsen_gene_corr_df[order(drugsen_gene_corr_df$r),]
    
    return(drugsen_gene_corr_df)

}


discovery.factor.omic.cor_plot <- function(omic_df, factor_df, factor_sel, type, gene, m='spearman',log=FALSE,just_plot_df=FALSE){



    if(type=='MOFA'){
        df_for_cal = merge((factor_df%>%filter(factor==factor_sel))[,c('sample','value')], omic_df)
    }else if(type=='PCA'){
        df_for_cal = merge(factor_df[,c('sample',factor_sel)],omic_df)
        colnames(df_for_cal)[2] = 'value'
    }else{
        print('ERROR type')
    }
    if(just_plot_df){
        return(df_for_cal)
    }
    

    p = ggpubr::ggscatter(df_for_cal, x=gene, y='value',conf.int = T, cor.method = m, add ='reg.line', cor.coef = T) + ylab(factor_sel)

    return(p)

}

# 生成带有组学mat的结果表
# df_for_cal由discovery.factor.omic.cor_plot生成just_plot_df=TRUE
# drugsen_gene_corr_df是discovery.factor.omic.cor 的结果表
discovery.factor.omic.cor.get_res_with_count_mat <- function(df_for_cal,drugsen_gene_corr_df){
    # ===== add count

    # df_for_cal = all_res_ls[['367C']]$phosph$df_for_cal
    rownames(df_for_cal) = NULL
    drugsen_gene_corr_df_with_count = drugsen_gene_corr_df
    drugsen_gene_corr_df_with_count = merge(drugsen_gene_corr_df_with_count,
                            data.frame(count=apply(df_for_cal[,c(3:ncol(df_for_cal))],2,function(x) sum(!is.na(x))))%>%rownames_to_column('gene'),
                            by='gene')

    drugsen_gene_corr_df_with_count = drugsen_gene_corr_df_with_count[order(drugsen_gene_corr_df_with_count$p,decreasing=FALSE),]



    # ===== add group
    drugsen_gene_corr_df_with_count_with_group = rbind(drugsen_gene_corr_df_with_count['NU',],
                                drugsen_gene_corr_df_with_count)
    drugsen_gene_corr_df_with_count_with_group[1,] = 'value'

    drugsen_gene_corr_df_with_count_with_group$gene = factor(drugsen_gene_corr_df_with_count_with_group$gene,
                                    level=drugsen_gene_corr_df_with_count_with_group$gene)



    # ===== add sample matrix
    drugsen_gene_corr_df_with_count_with_group = merge(
        drugsen_gene_corr_df_with_count_with_group,
        ( data.frame(t(df_for_cal%>%column_to_rownames('sample')))%>%rownames_to_column('gene') )[,],
        by='gene'
    )
    drugsen_gene_corr_df_with_count_with_group = 
        drugsen_gene_corr_df_with_count_with_group[order(drugsen_gene_corr_df_with_count_with_group$gene),]

   return(list(drugsen_gene_corr_df_with_count_with_group=drugsen_gene_corr_df_with_count_with_group,
              drugsen_gene_corr_df_with_count=drugsen_gene_corr_df_with_count))
}

############
## Pathway

discovery.factor.gsea <- function(cor_res, gmtfile, pvalueCutoff = 0.05){
    enrich_res = list()
    for(gg in names(gmtfile)){
        enrich_res[[gg]] = GSEA_enrich( (cor_res)[,c('gene','r')], gmtfile[[gg]] , pvalueCutoff = pvalueCutoff )
    }
    return(enrich_res)
}

############
## Volcano
# 可以自动挑选top来标注  《-------------------------
discovery.factor.volcano <- function(cor_res, show_genes=NULL, thr=0.25, thp=0.05, adj='p.adj',label_col='gene'){

    df_vol = cor_res

    colnames(df_vol)[1] = 'gene'

    df_vol$change = 'no'
    try_res = try({df_vol[(df_vol$r< -thr)&(df_vol[[adj]]< thp),]$change = 'neg'})
    try_res = try({df_vol[(df_vol$r> thr)&(df_vol[[adj]]< thp),]$change = 'pos'})

    df_vol$p_show = df_vol[[adj]]

    # ===
    p <- ggplot(
        # 数据、映射、颜色
        df_vol, aes(x = r, y = -log10(p_show), colour=change)) +
        geom_point(alpha=0.4, size=1) +
        scale_color_manual(values=c('blue',"grey","#ff4757"))+
        # 辅助线
        geom_vline(xintercept=c(-thr, thr),lty=4,col="black",lwd=0.8) +
        geom_hline(yintercept = -log10(thp),lty=4,col="black",lwd=0.8) +
        # 坐标轴
        labs(x="r",
            y=sprintf("-log10(%s)",adj))+
        theme_bw()+
        # 图例
        theme(plot.title = element_text(hjust = 0.5), 
                legend.position="right", 
                legend.title = element_blank())

    if(!is.null(show_genes)){
        p = p + ggrepel::geom_label_repel(data = df_vol%>%filter(gene%in%show_genes), aes(x = r, y = -log10(p_show), 
                        label = (df_vol%>%filter(gene%in%show_genes))$gene ),
                        size = 5, box.padding = unit(0.5, "lines"),
                        point.padding = unit(0.8, "lines"), 
                        segment.color = "black", 
                        show.legend = FALSE) +
            geom_point(data = df_vol%>%filter(gene%in%show_genes),
                    aes(x = r, y = -log10(p_show), label = gene),
                    alpha=0.4, size=3, color='black') +
            theme (legend.position = "none")
    }

    return(p)

}


############
## drug

discovery.factor.drug.cor <- function(drug_sen_df, factor_df, factor_sel, type, m='spearman'){

    if(type=='MOFA'){
        factor_df = reshape2::dcast(factors, sample~factor, value.var='value')
        if(factor_sel=='all'){
            factor_sel = colnames(factor_df)[-1]
        }
        df_for_cal = merge(factor_df[,c('sample',factor_sel)], drug_sen_df)
    }else if(type=='PCA'){
        if(factor_sel=='all'){
            factor_sel = colnames(factor_df)[-1]
        }
        df_for_cal = merge(factor_df[,c('sample',factor_sel)],drug_sen_df)
        # colnames(df_for_cal)[2] = 'value'
    }else{
        print('ERROR type')
    }
    

    cor_drug_factor_df_all_f = data.frame()

    for(f in factor_sel){
        # drug_cor_df = merge(factors,drug_sen_df)%>%filter(factor==factor_tmp)
        cor_drug_factor_df = list(factor=c(),inhibitor=c(),r=c(),p=c())

        for(drug_tmp in unique(df_for_cal$inhibitor)){
            # View(df_for_cal)
            df_cal = df_for_cal%>%filter(inhibitor==drug_tmp)
            df_cal = na.omit(df_cal)
            if(dim(df_cal)[1]<3) next
            cor_res= cor.test(as.formula(sprintf('~%s+ic50',f)),data=df_cal,method='spearman')

            cor_drug_factor_df$factor = c(cor_drug_factor_df$factor, f)
            cor_drug_factor_df$inhibitor = c(cor_drug_factor_df$inhibitor, drug_tmp)
            cor_drug_factor_df$r = c(cor_drug_factor_df$r, cor_res$estimate[[1]])
            cor_drug_factor_df$p = c(cor_drug_factor_df$p, cor_res$p.value)

        }

        cor_drug_factor_df = data.frame(cor_drug_factor_df)
        cor_drug_factor_df$adj.p = p.adjust(cor_drug_factor_df$p)
        cor_drug_factor_df = cor_drug_factor_df[order(cor_drug_factor_df$p),]

        cor_drug_factor_df_all_f = rbind(cor_drug_factor_df_all_f, cor_drug_factor_df)

    }

    return(cor_drug_factor_df_all_f)

}


discovery.factor.drug.plot <- function(drug_sen_df, factor_df, factor_sel, type, drug, m='spearman', log=TRUE){

    if(type=='MOFA'){
        df_for_cal = merge((factor_df%>%filter(factor==factor_sel))[,c('sample','value')], drug_sen_df)
    }else if(type=='PCA'){
        df_for_cal = merge(factor_df[,c('sample',factor_sel)],drug_sen_df)
        colnames(df_for_cal)[2] = 'value'
    }else{
        print('ERROR type')
    }

    if(log){
        df_for_cal[['ic50']] = log2(df_for_cal[['ic50']])

        p = ggpubr::ggscatter(df_for_cal%>%filter(inhibitor==drug), 
                x='value', y='ic50',add='reg.line',cor.coef = T, cor.method = m) + 
                xlab(factor_sel)+ylab(sprintf('log2(ic50)'))+labs(title=drug)
    }else{

        p = ggpubr::ggscatter(df_for_cal%>%filter(inhibitor==drug), 
                x='value', y='ic50',add='reg.line',cor.coef = T, cor.method = m) + 
                xlab(factor_sel)+ylab(sprintf('ic50'))+labs(title=drug)
    }

    return(p)
}

############
## mut

discovery.factor.mut.comp <- function(mut_df, factor_df, factor_sel, type, m='spearman'){

    if(type=='MOFA'){
        df_for_cal = merge((factor_df%>%filter(factor==factor_sel))[,c('sample','value')], mut_df)
    }else if(type=='PCA'){
        df_for_cal = merge(factor_df[,c('sample',factor_sel)],mut_df)
        colnames(df_for_cal)[2] = 'value'
    }else{
        print('ERROR type')
    }

    #df_cor_cal = merge(factor_df[,c('sample','Dim.2')],mut_df)
    #colnames(df_cor_cal)[2] = 'value'


    df_res = list(gene=c(),p=c(),fc=c())

    for(gene_tmp in colnames(df_cor_cal)[3:ncol(df_cor_cal)] ){
        val_mut = df_cor_cal[df_cor_cal[[gene_tmp]]==1,]$value
        val_wt = df_cor_cal[df_cor_cal[[gene_tmp]]==0,]$value
        if(length(val_mut)<3) next
        df_res$gene = c(df_res$gene,gene_tmp)
        df_res$p = c(df_res$p,wilcox.test(val_mut,val_wt)$p.val)
        df_res$fc = c(df_res$fc,median(val_mut)/median(val_wt))
    }
    df_res = data.frame(df_res)

    return(df_res)

}


discovery.factor.mut.plot <- function(mut_df, factor_df, factor_sel, type, gene, m='spearman'){
    if(type=='MOFA'){
        df_for_cal = merge((factor_df%>%filter(factor==factor_sel))[,c('sample','value')], mut_df)
    }else if(type=='PCA'){
        df_for_cal = merge(factor_df[,c('sample',factor_sel)],mut_df)
        colnames(df_for_cal)[2] = 'value'
    }else{
        print('ERROR type')
    }

    p=ggpubr::ggboxplot(df_for_cal,x=gene,y='value',add='jitter',color=gene)+stat_compare_means()

    return(p)

}

############
## survival

discovery.factor.survial.plot <- function(surv_df, factor_df, factor_sel, type){

    if(type=='MOFA'){
        df_for_cal = merge((factor_df%>%filter(factor==factor_sel))[,c('sample','value')], surv_df)
    }else if(type=='PCA'){
        df_for_cal = merge(factor_df[,c('sample',factor_sel)],surv_df)
        colnames(df_for_cal)[2] = 'value'
    }else{
        print('ERROR type')
    }

    p=plot_survival_contin(df_for_cal, sprintf('%s','value'), plot=T, title=NULL, OS_RF='OS')$plot+labs(title=sprintf('%s',factor_sel))
    
    return(p)

}

# factor_sel: can be alll
discovery.factor.survial.cal <- function(surv_df, factor_df, factor_sel='all', type){



    if(type=='MOFA'){
        factor_df = reshape2::dcast(factors, sample~factor, value.var='value')
        if(factor_sel=='all'){
            factor_sel = colnames(factor_df)[-1]
        }
        df_for_cal = merge(factor_df[,c('sample',factor_sel)], surv_df)
    }else if(type=='PCA'){
        if(factor_sel=='all'){
            factor_sel = colnames(factor_df)[-1]
        }
        df_for_cal = merge(factor_df[,c('sample',factor_sel)],surv_df)
        # colnames(df_for_cal)[2] = 'value'
    }else{
        print('ERROR type')
    }

    coxph_res_df = data.frame()
    for(f in factor_sel){
        coxph_res_df =  rbind(
            coxph_res_df,
            utils.do_coxph( df_for_cal, var=sprintf('%s',f),OS_RF='OS' )$coxph_df
        )
        
    }

    km_res_df = data.frame()
    for(f in factor_sel){
        km_res_df =  rbind(
            km_res_df,
            utils.do_km_survival( df_for_cal, var=sprintf('%s',f),OS_RF='OS' )$km_df
        )
        
    }

    return( list(coxph_res_df=coxph_res_df,km_res_df=km_res_df) )

}

# 独立预后因素




############
## Heatmap
# 要包括突变，生存，表达的热图，用comxplexheatmap


discovery.factor.heatmap <- function(omic_df, factor_df, factor_sel, gene_sel=NULL){

    if(is.matrix(omic_df)){
        omic_df = data.frame(omic_df)%>%rownames_to_column('sample')
    }

    tmp_df = factor_df[,c('sample',factor_sel)]
    tmp_df = tmp_df[order(tmp_df[[factor_sel]]),]


    omic_mat = omic_df%>%column_to_rownames('sample')
    if(!is.null(gene_sel)){
        omic_mat = omic_mat[,intersect(colnames(omic_mat),gene_sel)]
    }

    omic_mat = omic_mat[tmp_df$sample,]

    p = pheatmap(omic_mat,cluster_rows = F,border_color = NA )

    return(p)

}











############################################
## PREDICT ANALYSIS
############################################ 


predict.survival <- function(cor_res, gmtfile){



}




