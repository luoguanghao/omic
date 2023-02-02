# 这是一个用来测试feature的一个流程，这个流程
## ~计算ssGSEA分数
## ~评估ssGSEA分数指标
## ~对比别的指标：BLAST，NPM1，单基因指标:输出表格
## 
## 
##
## ~~~需要进行的调整：roc计算免输出，ROC的参数调整
## compare函数需要加入计算相关系数的模块

library(stringr)
library(ggpubr) 
library(ggplot2)
library(GSEABase)
#library(GSVAdata)
#data(c2BroadSets)
#c2BroadSets
 
library(Biobase)
#library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)
library(tidyverse)
library(pheatmap)
#library(DESeq2)

# library(edgeR)
library(limma)
library(pROC)

library(patchwork)


source('/home/lgh/my_project/omic/modules/plot/boxplot.R')
# source('/mnt/d/omic/modules/plot/auc.R')
source('/home/lgh/my_project/omic/modules/plot/survival.R')

## get and show overlap
# some gene_ls combined in list format
# 
# input: all_gene_list, num_overlap
# 获取一个矩阵，每列是一个药，每行是一个基因，描述了每个药都有哪些基因被选中
# output: overlap_gene, maybe 3 overlap or 4 overlap
# 
get_plot_overlap <- function(all_gene_list, num_overlap){

    all_gene = c()
    for(i in names(all_gene_list)){
        all_gene = c(all_gene, all_gene_list[[i]])
    }
    all_gene = unique(all_gene)

    drug_ls = c()
    value_ls = c()
    for(i in names(all_gene_list)){
        drug_ls = c( drug_ls, rep(i,length(all_gene)) )
        value_ls = c( value_ls, all_gene%in%all_gene_list[[i]] )
    }

    df_for_plot = as.data.frame(list( gene=rep(all_gene,4),
                                    drug=drug_ls,
                                    value=value_ls
                                    ))

    df_mean=group_by(df_for_plot%>%filter(value==TRUE), gene) %>% summarise(count = n())
    df_for_plot$gene = factor(df_for_plot$gene, levels = as.vector(df_mean[order(df_mean$count),][['gene']]))

    cols=c('TRUE'='red','FALSE'='black')

    p = df_for_plot%>%ggplot(aes(x=drug,y=gene))+
    geom_tile(aes(fill=value),color="white",size=1)+ #color和size分别指定方块边线的颜色和粗细
    scale_x_discrete("",expand = c(0,0))+ #不显示横纵轴的label文本；画板不延长
    scale_y_discrete("",expand = c(0,0))+
    coord_equal(0.1) +
    scale_fill_manual(values = cols)+ #指定自定义的颜色
    theme(
        axis.text.x.bottom = element_text(size=10),axis.text.y.left = element_text(size = 12), #修改坐标轴文本大小
        axis.ticks = element_blank(), #不显示坐标轴刻度
        legend.title = element_blank() #不显示图例title
    )
    plot(p)
    
    overlap_gene = as.vector( (df_mean%>%filter(count>=num_overlap))$gene )

    return(overlap_gene)

}


# get_sen_unsen
## list(sen_sample=sen_sample, unsen_sample=unsen_sample)
get_sen_unsen <- function(drug_sen_df, drug, sp_sample_ls, exp_mat,auc_ic50='auc',precent=0.3){

    # gene_ls = intersect(gene_ls,rownames(exp_mat))
    ## divide
    ### get sen unsen

    sp_drug_sen_df = drug_sen_df[drug_sen_df['inhibitor']==drug,]%>%filter(lab_id%in%sp_sample_ls)
    sp_drug_sen_df = sp_drug_sen_df[order(sp_drug_sen_df[[auc_ic50]]),]
    sp_drug_sen_df$lab_id = factor(sp_drug_sen_df$lab_id,level=as.vector(sp_drug_sen_df$lab_id))

    top_bottom_number = round(dim(sp_drug_sen_df)[1]*precent)
    sen_sample = as.vector(sp_drug_sen_df[1:top_bottom_number,]$lab_id)
    unsen_sample = as.vector(sp_drug_sen_df[(dim(sp_drug_sen_df)[1]-top_bottom_number):dim(sp_drug_sen_df)[1],]$lab_id)
    
    return( list(sen_sample=sen_sample, unsen_sample=unsen_sample) )
}


# ssgsea_score_gene_set
# exp_mat:表达矩阵，df格式
# gene_ls需要计算ssgsea分数的vector
ssgsea_score_gene_set <- function(exp_mat, gene_ls){
    gene_ls = intersect(gene_ls,rownames(exp_mat))

    gene_list = list(gene_ls=gene_ls)
    gsva_matrix <- gsva(as.matrix(exp_mat[,]), 
                        gene_list,
                        method='ssgsea', 
                        kcdf='Gaussian', 
                        abs.ranking=TRUE, 
                        parallel.sz=1)
    #return(gsva_matrix)
    df_ssgsea_score = data.frame(sample=colnames(gsva_matrix), 
                            ssGSEA_score=as.vector(gsva_matrix[1,]))

    return(df_ssgsea_score)
}





# test_gene_ls_in_drug() 测试feature gene的流程
#
# 分析基因list的ssgsea分数与药敏的关系
## pheatmap,boxplot,roc,scatter
##############
#
# input:
## drug_sen_df:inhibitor;lab_id;ic50;auc
## drug
## sp_sample_ls:有用的样本
## exp_mat：表达矩阵：df:rowname is gene, colnames is sample
## gene_ls:挑选出来拿来ssGSEA的基因
## sen_unsen_list: $sen_sample $unsen_sample
##  auc_ic50='auc' 'ic50'
# 带有置信区间
test_gene_ls_in_drug <- function(drug_sen_df, drug, sp_sample_ls, sen_unsen_list, exp_mat, gene_ls, 
                                        auc_ic50='auc', plot_flag=FALSE){
    # drug_sen_df
    # sp_sample_ls
    # exp_mat
    # gene_ls
    #
    gene_ls = intersect(gene_ls,rownames(exp_mat))
    ## divide
    ### get sen unsen

    sp_drug_sen_df = drug_sen_df[drug_sen_df['inhibitor']==drug,]%>%filter(lab_id%in%sp_sample_ls)
    sp_drug_sen_df = sp_drug_sen_df[order(sp_drug_sen_df[[auc_ic50]]),]
    sp_drug_sen_df$lab_id = factor(sp_drug_sen_df$lab_id,level=as.vector(sp_drug_sen_df$lab_id))

    sen_sample = sen_unsen_list$sen_sample
    unsen_sample = sen_unsen_list$unsen_sample

    #top_bottom_number = round(dim(sp_drug_sen_df)[1]*0.3)
    #sen_sample = as.vector(sp_drug_sen_df[1:top_bottom_number,]$lab_id)
    #unsen_sample = as.vector(sp_drug_sen_df[(dim(sp_drug_sen_df)[1]-top_bottom_number):dim(sp_drug_sen_df)[1],]$lab_id)

    # par(mfrow=c(2,2))
    ## 看看overlap基因的热图
    ### pheatmap
    list_for_pheatmap = list(
        exp_mat=exp_mat[gene_ls, c(sen_sample,unsen_sample)],
        anno_nc_df=data.frame(type=c(rep('sen',length(sen_sample)),rep('unsen',length(unsen_sample))),row.names=c(sen_sample,unsen_sample)),
        drug=drug
    )
    if(plot_flag==TRUE){
        anno_nc_df = data.frame(type=c(rep('sen',length(sen_sample)),rep('unsen',length(unsen_sample))),row.names=c(sen_sample,unsen_sample))
        phm = pheatmap(
            exp_mat[gene_ls, c(sen_sample,unsen_sample)],
            cluster_cols = FALSE,
            #cluster_rows = FALSE,
            scale = 'row',
            color = colorRampPalette(c("green","black", "red"))(50),
            # filename='./cencluster_hmap1.png',
            cellwidth = 8, cellheight = 18,border=FALSE,
            #cutree_cols=2, 
            annotation_col  = anno_nc_df,
            main = drug
            )
    }
    ## ssGSEA打分

    gene_list = list(gene_ls=gene_ls)
    gsva_matrix <- gsva(as.matrix(exp_mat[,]), 
                        gene_list,
                        method='ssgsea', 
                        kcdf='Gaussian', 
                        abs.ranking=TRUE, 
                        parallel.sz=1)
    #return(gsva_matrix)
    df_ssgsea_score = data.frame(sample=colnames(gsva_matrix), ssGSEA_score=as.vector(gsva_matrix[1,]))
    df_plot = as.data.frame(gsva_matrix[,c(sen_sample,unsen_sample)])
    df_plot$sample = rownames(df_plot)
    df_plot$sen_unsen = c(rep('sen',length(sen_sample)), rep('unsen',length(unsen_sample)))

    colnames(df_plot) = c('ssGSEA_score','sample','sen_unsen')
    
        # plot
    list_for_boxplot_roc_scatter = list(
        df_plot=df_plot
    )
    if(plot_flag==TRUE){
        #p=ggplot(data=df_plot,aes(x=sen_unsen, y=gene_ls)) + geom_boxplot(aes(fill=sen_unsen)) + 
        #    stat_compare_means() + 
        #    labs(title=drug) + ylab('ssGSEA score overlap ND gene')+
        #    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
        #plot(p)
        p_boxplot = my_boxplot(df_plot, x='sen_unsen', y='ssGSEA_score', 
                                label_list=list(labs=drug,ylab='ssGSEA Score',xlab=''), 
                                para_list=NA, out_name=NA)
    }
    # roc(response, predictor)
    roc1 <- roc(df_plot$sen_unsen, df_plot$ssGSEA_score,ci=TRUE, boot.n=10000)

    list_for_boxplot_roc_scatter$roc=roc1
    list_for_boxplot_roc_scatter$df_scatter=merge(df_ssgsea_score, sp_drug_sen_df, by.x='sample', by.y='lab_id',drug=drug)
    if(plot_flag==TRUE){
        plot(roc1,
            thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
            print.thres="best",
            print.auc=TRUE)
        # plot(proc)

        # colnames(df_plot) = c('ssGSEA_score','sample','sen_unsen')

        pscatter=ggscatter(merge(df_ssgsea_score, sp_drug_sen_df, by.x='sample', by.y='lab_id'), 
                    x=auc_ic50, y='ssGSEA_score', palette = "nature", 
                    add = "reg.line", conf.int = TRUE,title=drug)+stat_cor()

        print(pscatter)
    }else{
        p_boxplot=NULL
        pscatter=NULL
    }

    return(
        list(gsva_matrix=gsva_matrix,
            df_ssgsea_score=df_ssgsea_score,
            df_plot=df_plot,
            df_for_contin_plot=merge(df_ssgsea_score, sp_drug_sen_df, by.x='sample', by.y='lab_id'),
            sp_drug_sen_df=sp_drug_sen_df,
            p_bp=p_boxplot,
            pscatter=pscatter,
            roc=roc1,
            list_for_boxplot_roc_scatter=list_for_boxplot_roc_scatter,
            list_for_pheatmap=list_for_pheatmap)
        )

}



####
# 只做ssGSEA的分析，计算每个样本的SSGSEA分数，并且标上sen，unsen的标签，与上一个函数相比不做图
# 带有置信区间
test_gene_ls_in_drug_for_Exhaustion <- function(drug_sen_df, drug, sp_sample_ls, sen_unsen_list, exp_mat, gene_ls, 
                                                auc_ic50='auc',verbose=FALSE,parallel.sz=1){
    # drug_sen_df
    # sp_sample_ls
    # exp_mat
    # gene_ls
    #
    gene_ls = intersect(gene_ls,rownames(exp_mat))
    ## divide
    ### get sen unsen


    sen_sample = sen_unsen_list$sen_sample
    unsen_sample = sen_unsen_list$unsen_sample

    #top_bottom_number = round(dim(sp_drug_sen_df)[1]*0.3)
    #sen_sample = as.vector(sp_drug_sen_df[1:top_bottom_number,]$lab_id)
    #unsen_sample = as.vector(sp_drug_sen_df[(dim(sp_drug_sen_df)[1]-top_bottom_number):dim(sp_drug_sen_df)[1],]$lab_id)

    ## ssGSEA打分

    gene_list = list(gene_ls=gene_ls)
    gsva_matrix <- gsva(as.matrix(exp_mat[,]), 
                        gene_list,
                        method='ssgsea', 
                        kcdf='Gaussian', 
                        abs.ranking=TRUE, 
                        parallel.sz=parallel.sz,
                        verbose=verbose)
    #return(gsva_matrix)
    df_ssgsea_score = data.frame(sample=colnames(gsva_matrix), ssGSEA_score=as.vector(gsva_matrix[1,]))
    df_plot = as.data.frame(gsva_matrix[,c(sen_sample,unsen_sample)])
    df_plot$sample = rownames(df_plot)
    df_plot$sen_unsen = c(rep('sen',length(sen_sample)), rep('unsen',length(unsen_sample)))

    colnames(df_plot) = c('ssGSEA_score','sample','sen_unsen')

    return(
        list(gsva_matrix=gsva_matrix,
            df_ssgsea_score=df_ssgsea_score,
            df_plot=df_plot)
        )
}




#
#
# drug_sen_df: inhibitor;lab_id;ic50;auc
# blastbm_df: lab_id;Blasts.in.BM
# 带有置信区间
compare_sen_unsen_blast <- function(drug_sen_df, drug, sp_sample_ls, sen_unsen_list, blastbm_df, aucic50, pb_bm, plot_flag=FALSE){

    #gene_ls = intersect(gene_ls,rownames(exp_mat))
    ## divide
    ### get sen unsen

    sp_drug_sen_df = drug_sen_df[drug_sen_df['inhibitor']==drug,]%>%filter(lab_id%in%sp_sample_ls)
    sp_drug_sen_df = sp_drug_sen_df[order(sp_drug_sen_df[[aucic50]]),]
    sp_drug_sen_df$lab_id = factor(sp_drug_sen_df$lab_id,level=as.vector(sp_drug_sen_df$lab_id))

    #top_bottom_number = round(dim(sp_drug_sen_df)[1]*0.3)
    #sen_sample = as.vector(sp_drug_sen_df[1:top_bottom_number,]$lab_id)
    #unsen_sample = as.vector(sp_drug_sen_df[(dim(sp_drug_sen_df)[1]-top_bottom_number):dim(sp_drug_sen_df)[1],]$lab_id)

    sen_sample = sen_unsen_list$sen_sample
    unsen_sample = sen_unsen_list$unsen_sample

    # mark sen unsen
    sp_drug_sen_df$sen_unsen = NA
    sp_drug_sen_df[sp_drug_sen_df$lab_id%in%sen_sample,]$sen_unsen='sen'
    sp_drug_sen_df[sp_drug_sen_df$lab_id%in%unsen_sample,]$sen_unsen='unsen'
    colnames(sp_drug_sen_df)[length(colnames(sp_drug_sen_df))] = 'sen_unsen'


    # prepare df_for_plot
    # drug_sen_flag; blastbm; survival; sample
    df_for_plot = na.omit(merge(blastbm_df,sp_drug_sen_df,by='lab_id'))
    if(pb_bm=='Blasts.in.BM'){
        df_for_plot = df_for_plot%>%filter(is.na(sen_unsen)==FALSE)%>%filter(Blasts.in.BM!='')
    }else{
        df_for_plot = df_for_plot%>%filter(is.na(sen_unsen)==FALSE)%>%filter(Blasts.in.PB!='')
    }
    df_for_plot = df_for_plot[,c('lab_id',pb_bm,'sen_unsen')]
    colnames(df_for_plot) = c('sample','y_val','x_val')
    if(plot_flag==TRUE){
        my_boxplot(df_for_plot, label_list=list(labs=drug,ylab=pb_bm,xlab='sen_unsen'))
    }  
    roc1 <- roc(df_for_plot$x_val, df_for_plot$y_val, 
                ci=TRUE, boot.n=10000) # 先senunsen,再blast
    if(plot_flag==TRUE){
        plot(roc1,
            thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
            print.thres="best",
            print.auc=TRUE)

        p = ggscatter(merge(blastbm_df,drug_sen_df%>%filter(inhibitor==drug)),
                x=pb_bm,y=aucic50,
                add = "reg.line",conf.int = TRUE,
                title=drug,ylab='IC50' ,xlab=pb_bm )+stat_cor()
        plot(p)
    }

    #pscatter=ggscatter(merge(df_ssgsea_score, sp_drug_sen_df, by.x='sample', by.y='lab_id'), 
    #            x=auc_ic50, y='ssGSEA_score', palette = "nature", 
    #            add = "reg.line", conf.int = TRUE,title=drug)+stat_cor()
    if(pb_bm=='Blasts.in.BM'){
        colnames(df_for_plot) = c('sample','%Blasts.in.BM','sen_unsen')
    }else{
        colnames(df_for_plot) = c('sample','%Blasts.in.PB','sen_unsen')
    }
    
    return(list(
        df_for_binary_plot = df_for_plot,
        df_for_contin_plot = merge(blastbm_df,drug_sen_df%>%filter(inhibitor==drug)),
        roc = roc1
    ))
}

# 
## drugsen_unsen_list: $sen_sample, $unsen_sample
## marker_df: $lab_id, $marker
## 带有置信区间
compare_binary_marker <- function(drugsen_unsen_list,marker_df,plot_flag=FALSE){

    df_for_binary_marker = data.frame(lab_id=c(drugsen_unsen_list$sen_sample,drugsen_unsen_list$unsen_sample),
                                    sen_unsen=c( rep('sen',length(drugsen_unsen_list$sen_sample)),rep('unsen',length(drugsen_unsen_list$unsen_sample)) )
                                    )
    df_for_binary_marker = merge(df_for_binary_marker,marker_df, by='lab_id')

    df_for_binary_marker$marker[df_for_binary_marker$marker=='wt'] = 0
    df_for_binary_marker$marker[df_for_binary_marker$marker=='mut'] = 1
    df_for_binary_marker$marker = as.numeric(df_for_binary_marker$marker)


    roc1 <- roc(df_for_binary_marker$sen_unsen, 
                df_for_binary_marker$marker,ci=TRUE,boot.n=10000)
    if(plot_flag==TRUE){
        plot(roc1,
            thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
            print.thres="best",
            print.auc=TRUE)
    }
    return(list(df_for_binary_marker=df_for_binary_marker,roc=roc1))

}


#### 
# dir,filename,drug
# exp_mat 表达矩阵 gene*sample, 
# sp_sample_ls:特定的sample的vector，maybe 'itd'
# auc_ic50
# ~~~sen_unsen_list
# drug_sen_df: colnames = c('inhibitor','lab_id','ic50','auc')
#
# blastbm_df: colnames = c('lab_id','Blasts.in.BM')
# NPM1_df: colnames = c('sample','NPM1_status')
# ~~~hl_8gene_set
# cebpa_set: 感兴趣的基因list
#
# 需要最后给出一个表格：带有置信区间
compare_3_metric <- function(dir,filename,drug,exp_mat,sp_sample_ls,
                    drug_sen_df,blastbm_df,NPM1_df,cebpa_set,
                    senunsen_prec=0.3,auc_ic50='ic50',plot_flag=FALSE) {

    # genels_name = 'newcebpaMQC'
    # result_df = list(hl8=c(),cebpa=c(),blast=c(),npm1=c())

    #  drug = our_drugs[1] # 
    print(drug)
    ## choose sp sample
    # itd_exp_drug_sample = intersect(intersect(ITD_sample_ls,colnames(exp_mat)), (drug_sen_df%>%filter(inhibitor==drug))$lab_id)

    # list(sen_sample=sen_sample, unsen_sample=unsen_sample)
    get_sen_unsen_res = get_sen_unsen(drug_sen_df=drug_sen_df, drug=drug,
                            sp_sample_ls=itd_exp_drug_sample, exp_mat=exp_mat, 
                            auc_ic50=auc_ic50,precent=senunsen_prec) ##

    # test for ND overlap gene
    #res_df_hl8 = test_gene_ls_in_drug(drug_sen_df=drug_sen_df, drug=drug,
    #                    sp_sample_ls=sp_sample_ls, sen_unsen_list=get_sen_unsen_res, #
    #                    exp_mat=exp_mat, gene_ls=hl_8gene_set,auc_ic50='ic50')

    res_df_cebpa = test_gene_ls_in_drug(drug_sen_df=drug_sen_df, drug=drug,
                        sp_sample_ls=sp_sample_ls, sen_unsen_list=get_sen_unsen_res, #
                        exp_mat=exp_mat, gene_ls=cebpa_set,auc_ic50=auc_ic50)

    blast_res_list = compare_sen_unsen_blast(drug_sen_df=drug_sen_df, 
                                            drug=drug, 
                                            sp_sample_ls=sp_sample_ls,
                                            sen_unsen_list=get_sen_unsen_res, #
                                            blastbm_df=blastbm_df,
                                            aucic50=auc_ic50,pb_bm='Blasts.in.BM')

    marker_df = NPM1_df
    colnames(marker_df) = c('lab_id','marker')
    marker_df$marker = as.character(marker_df$marker)
    npm1_res_list = compare_binary_marker(drugsen_unsen_list=get_sen_unsen_res,
                        marker_df=marker_df)

    #roc1_2 <- roc(res_df_cebpa$df_plot$sen_unsen, 
    #            res_df_cebpa$df_plot$ssGSEA_score,ci=TRUE,boot.n=10000)
    roc1_2 = res_df_cebpa$roc
    #roc2 <- roc(blast_res_list$df_for_binary_plot$x_val, 
    #            blast_res_list$df_for_binary_plot$y_val,ci=TRUE,boot.n=10000)
    roc2 = blast_res_list$roc
    #roc3 <- roc(npm1_res_list$df_for_binary_marker$sen_unsen, 
    #            npm1_res_list$df_for_binary_marker$marker,ci=TRUE,boot.n=10000)
    roc3 = npm1_res_list$roc
    #cat('hl8',roc1$auc,'\n')
    #cat('cebpa',roc1_2$auc,'\n')
    #cat('blast',roc2$auc,'\n')
    #cat('npm1',roc3$auc,'\n')
    #result_df$hl8 = c(result_df$hl8,roc1$auc)
    result_df = list(drug=c(),cebpa=c(),cebpa_ci=c(),blast=c(),blast_ci=c(),
                        cebpa_blast_p=c(),npm1=c(),npm1_ci=c(),cebpa_npm1_p=c())
    
    result_df$drug = c(result_df$drug,drug)
    
    result_df$cebpa = c(result_df$cebpa,roc1_2$auc)
    result_df$blast = c(result_df$blast,roc2$auc)
    result_df$npm1 = c(result_df$npm1,roc3$auc)
    
    result_df$cebpa_blast_p = c(result_df$cebpa_blast_p,roc.test(roc1_2,roc2, reuse.auc=TRUE,boot.n=10000)$p.value)
    result_df$cebpa_npm1_p = c(result_df$cebpa_npm1_p,roc.test(roc1_2,roc3, reuse.auc=TRUE,boot.n=10000)$p.value)
    
    result_df$cebpa_ci = c(result_df$cebpa_ci,paste(sprintf('%.2f',roc1_2$ci),collapse=('_')))
    result_df$blast_ci = c(result_df$blast_ci,paste(sprintf('%.2f',roc2$ci),collapse=('_')))
    result_df$npm1_ci = c(result_df$npm1_ci,paste(sprintf('%.2f',roc3$ci),collapse=('_')))

    result_df = data.frame(result_df)

    if(plot_flag==TRUE){
        text = paste(
            c(paste('ssGSEA score=',sprintf('%.2f',roc1_2$auc),sep=''),
            paste('BLAST=',sprintf('%.2f',roc2$auc),sep=''),
            paste('NPM1=',sprintf('%.2f',roc3$auc),sep='')), collapse='\n'
                    )

        ggroc(list(ssGSEA_score_8_genes = roc1, ssGSEA_score_cebpa_genes = roc1_2, blast = roc2, npm1 = roc3))+labs(title=drug) + theme_bw()

        ggroc(list(ssGSEA_score_cebpa_genes = roc1_2, blast = roc2, npm1 = roc3))+ 
                annotate('text', x = 0.4, y = 0.10, alpha = .99, colour = "black",label = text,hjust=0) +
                theme_bw()+theme(text=element_text(size=15,face='bold',family='serif'))
        ggsave(sprintf('%s/%s_%s_%s.pdf', dir, filename, auc_ic50, drug), dpi=70, width=7.5, height=5) #
    }
    return(list(
        result_df=result_df,
        data=list(res_df_cebpa=res_df_cebpa,
                blast_res_list=blast_res_list,
                npm1_res_list=npm1_res_list)   
    ))

}




#### 
# dir,filename,drug
# exp_mat 表达矩阵 gene*sample, 
# sp_sample_ls:特定的sample的vector，maybe 'itd'
# auc_ic50
# sen_unsen_list: list(sen_sample=sen_sample, unsen_sample=unsen_sample)
# drug_sen_df: colnames = c('inhibitor','lab_id','ic50','auc')
#
# blastbm_df: colnames = c('lab_id','Blasts.in.BM')
# NPM1_df: colnames = c('sample','NPM1_status')
# ~~~hl_8gene_set
# cebpa_set: 感兴趣的基因list
#
# 需要最后给出一个表格：带有置信区间
#
#
# get_sen_unsen_res = get_sen_unsen(drug_sen_df=drug_sen_df, drug=drug,
#                        sp_sample_ls=itd_exp_drug_sample, exp_mat=exp_mat, 
#                        auc_ic50=auc_ic50,precent=senunsen_prec) ##
compare_3_metric_new <- function(dir,filename,drug,exp_mat,sp_sample_ls,
                    drug_sen_df,blastbm_df,NPM1_df,cebpa_set,
                    sen_unsen_list,auc_ic50='ic50',plot_flag=FALSE) {

    # genels_name = 'newcebpaMQC'
    # result_df = list(hl8=c(),cebpa=c(),blast=c(),npm1=c())

    #  drug = our_drugs[1] # 
    print(drug)
    ## choose sp sample
    # itd_exp_drug_sample = intersect(intersect(ITD_sample_ls,colnames(exp_mat)), (drug_sen_df%>%filter(inhibitor==drug))$lab_id)

    # list(sen_sample=sen_sample, unsen_sample=unsen_sample)
    #get_sen_unsen_res = get_sen_unsen(drug_sen_df=drug_sen_df, drug=drug,
    #                        sp_sample_ls=sp_sample_ls, exp_mat=exp_mat, 
    #                        auc_ic50=auc_ic50,precent=senunsen_prec) ##

    # test for ND overlap gene
    #res_df_hl8 = test_gene_ls_in_drug(drug_sen_df=drug_sen_df, drug=drug,
    #                    sp_sample_ls=sp_sample_ls, sen_unsen_list=get_sen_unsen_res, #
    #                    exp_mat=exp_mat, gene_ls=hl_8gene_set,auc_ic50='ic50')

    res_df_cebpa = test_gene_ls_in_drug(drug_sen_df=drug_sen_df, drug=drug,
                        sp_sample_ls=sp_sample_ls, sen_unsen_list=sen_unsen_list, #
                        exp_mat=exp_mat, gene_ls=cebpa_set,auc_ic50=auc_ic50)

    blast_res_list = compare_sen_unsen_blast(drug_sen_df=drug_sen_df, 
                                            drug=drug, 
                                            sp_sample_ls=sp_sample_ls,
                                            sen_unsen_list=sen_unsen_list, #
                                            blastbm_df=blastbm_df,
                                            aucic50=auc_ic50,pb_bm='Blasts.in.BM')

    marker_df = NPM1_df
    colnames(marker_df) = c('lab_id','marker')
    marker_df$marker = as.character(marker_df$marker)
    npm1_res_list = compare_binary_marker(drugsen_unsen_list=sen_unsen_list,
                        marker_df=marker_df)

    #roc1_2 <- roc(res_df_cebpa$df_plot$sen_unsen, 
    #            res_df_cebpa$df_plot$ssGSEA_score,ci=TRUE,boot.n=10000)
    roc1_2 = res_df_cebpa$roc
    #roc2 <- roc(blast_res_list$df_for_binary_plot$x_val, 
    #            blast_res_list$df_for_binary_plot$y_val,ci=TRUE,boot.n=10000)
    roc2 = blast_res_list$roc
    #roc3 <- roc(npm1_res_list$df_for_binary_marker$sen_unsen, 
    #            npm1_res_list$df_for_binary_marker$marker,ci=TRUE,boot.n=10000)
    roc3 = npm1_res_list$roc
    #cat('hl8',roc1$auc,'\n')
    #cat('cebpa',roc1_2$auc,'\n')
    #cat('blast',roc2$auc,'\n')
    #cat('npm1',roc3$auc,'\n')
    #result_df$hl8 = c(result_df$hl8,roc1$auc)
    result_df = list(drug=c(),cebpa=c(),cebpa_ci=c(),blast=c(),blast_ci=c(),
                        cebpa_blast_p=c(),npm1=c(),npm1_ci=c(),cebpa_npm1_p=c())
    
    result_df$drug = c(result_df$drug,drug)
    
    result_df$cebpa = c(result_df$cebpa,roc1_2$auc)
    result_df$blast = c(result_df$blast,roc2$auc)
    result_df$npm1 = c(result_df$npm1,roc3$auc)
    
    result_df$cebpa_blast_p = c(result_df$cebpa_blast_p,roc.test(roc1_2,roc2, reuse.auc=TRUE,boot.n=10000)$p.value)
    result_df$cebpa_npm1_p = c(result_df$cebpa_npm1_p,roc.test(roc1_2,roc3, reuse.auc=TRUE,boot.n=10000)$p.value)
    
    result_df$cebpa_ci = c(result_df$cebpa_ci,paste(sprintf('%.2f',roc1_2$ci),collapse=('_')))
    result_df$blast_ci = c(result_df$blast_ci,paste(sprintf('%.2f',roc2$ci),collapse=('_')))
    result_df$npm1_ci = c(result_df$npm1_ci,paste(sprintf('%.2f',roc3$ci),collapse=('_')))

    result_df = data.frame(result_df)

    if(plot_flag==TRUE){
        text = paste(
            c(paste('ssGSEA score=',sprintf('%.2f',roc1_2$auc),sep=''),
            paste('BLAST=',sprintf('%.2f',roc2$auc),sep=''),
            paste('NPM1=',sprintf('%.2f',roc3$auc),sep='')), collapse='\n'
                    )

        # ggroc(list(ssGSEA_score_8_genes = roc1, ssGSEA_score_cebpa_genes = roc1_2, blast = roc2, npm1 = roc3))+labs(title=drug) + theme_bw()

        p=ggroc(list(ssGSEA_score_cebpa_genes = roc1_2, blast = roc2, npm1 = roc3),size=1)+ 
                annotate('text', x = 0.4, y = 0.10, alpha = .99, colour = "black",label = text,hjust=0) +
                theme_bw()+theme(text=element_text(size=15,face='bold',family='serif'))
        ggsave(sprintf('%s/%s_%s_%s.pdf', dir, filename, auc_ic50, drug), dpi=70, width=7.5, height=5) #
    }
    return(list(
        result_df=result_df,
        data=list(res_df_cebpa=res_df_cebpa,
                blast_res_list=blast_res_list,
                npm1_res_list=npm1_res_list)   
    ))

}


###############################################################
# other comparation
###############################################################

compare_pb_bm_blast <- function(dir,filename,drug,exp_mat,sp_sample_ls,
                    drug_sen_df,blastbm_df,blastpb_df,NPM1_df,cebpa_set,
                    senunsen_prec=0.3,auc_ic50='ic50',plot_flag=FALSE){

    # genels_name = 'newcebpaMQC'
    # result_df = list(hl8=c(),cebpa=c(),blast=c(),npm1=c())

    #  drug = our_drugs[1] # 
    print(drug)
    ## choose sp sample
    # itd_exp_drug_sample = intersect(intersect(ITD_sample_ls,colnames(exp_mat)), (drug_sen_df%>%filter(inhibitor==drug))$lab_id)

    get_sen_unsen_res = get_sen_unsen(drug_sen_df=drug_sen_df, drug=drug,
                            sp_sample_ls=itd_exp_drug_sample, exp_mat=exp_mat, 
                            auc_ic50=auc_ic50,precent=senunsen_prec) ##

    # test for ND overlap gene
    #res_df_hl8 = test_gene_ls_in_drug(drug_sen_df=drug_sen_df, drug=drug,
    #                    sp_sample_ls=sp_sample_ls, sen_unsen_list=get_sen_unsen_res, #
    #                    exp_mat=exp_mat, gene_ls=hl_8gene_set,auc_ic50='ic50')

    res_df_cebpa = test_gene_ls_in_drug(drug_sen_df=drug_sen_df, drug=drug,
                        sp_sample_ls=sp_sample_ls, sen_unsen_list=get_sen_unsen_res, #
                        exp_mat=exp_mat, gene_ls=cebpa_set,auc_ic50=auc_ic50)

    blast_bm_res_list = compare_sen_unsen_blast(drug_sen_df=drug_sen_df, 
                                            drug=drug, 
                                            sp_sample_ls=sp_sample_ls,
                                            sen_unsen_list=get_sen_unsen_res, #
                                            blastbm_df=blastbm_df,
                                            aucic50=auc_ic50,pb_bm='Blasts.in.BM')

    blast_pb_res_list = compare_sen_unsen_blast(drug_sen_df=drug_sen_df, 
                                            drug=drug, 
                                            sp_sample_ls=sp_sample_ls,
                                            sen_unsen_list=get_sen_unsen_res, #
                                            blastbm_df=blastpb_df,
                                            aucic50=auc_ic50,pb_bm='Blasts.in.PB')

    marker_df = NPM1_df
    colnames(marker_df) = c('lab_id','marker')
    marker_df$marker = as.character(marker_df$marker)
    npm1_res_list = compare_binary_marker(drugsen_unsen_list=get_sen_unsen_res,
                        marker_df=marker_df)

    #roc1_2 <- roc(res_df_cebpa$df_plot$sen_unsen, 
    #            res_df_cebpa$df_plot$ssGSEA_score,ci=TRUE,boot.n=10000)
    roc1_2 = res_df_cebpa$roc
    #roc2 <- roc(blast_res_list$df_for_binary_plot$x_val, 
    #            blast_res_list$df_for_binary_plot$y_val,ci=TRUE,boot.n=10000)
    roc2_1 = blast_bm_res_list$roc
    roc2_2 = blast_pb_res_list$roc
    #roc3 <- roc(npm1_res_list$df_for_binary_marker$sen_unsen, 
    #            npm1_res_list$df_for_binary_marker$marker,ci=TRUE,boot.n=10000)
    roc3 = npm1_res_list$roc
    #cat('hl8',roc1$auc,'\n')
    #cat('cebpa',roc1_2$auc,'\n')
    #cat('blast',roc2$auc,'\n')
    #cat('npm1',roc3$auc,'\n')
    #result_df$hl8 = c(result_df$hl8,roc1$auc)
    result_df = list(drug=c(),cebpa=c(),cebpa_ci=c(),
                    blastbm=c(),blastbm_ci=c(),cebpa_blastbm_p=c(),
                    blastpb=c(),blastpb_ci=c(),cebpa_blastpb_p=c(),
                    npm1=c(),npm1_ci=c(),cebpa_npm1_p=c())
    
    result_df$drug = c(result_df$drug,drug)
    
    result_df$cebpa = c(result_df$cebpa,roc1_2$auc)
    result_df$blastbm = c(result_df$blastbm,roc2_1$auc)
    result_df$blastpb = c(result_df$blastpb,roc2_2$auc)
    result_df$npm1 = c(result_df$npm1,roc3$auc)
    
    result_df$cebpa_blastbm_p = c(result_df$cebpa_blastbm_p,roc.test(roc1_2,roc2_1, reuse.auc=TRUE,boot.n=10000)$p.value)
    result_df$cebpa_blastpb_p = c(result_df$cebpa_blastpb_p,roc.test(roc1_2,roc2_2, reuse.auc=TRUE,boot.n=10000)$p.value)
    result_df$cebpa_npm1_p = c(result_df$cebpa_npm1_p,roc.test(roc1_2,roc3, reuse.auc=TRUE,boot.n=10000)$p.value)
    
    result_df$cebpa_ci = c(result_df$cebpa_ci,paste(sprintf('%.2f',roc1_2$ci),collapse=('_')))
    result_df$blastbm_ci = c(result_df$blastbm_ci,paste(sprintf('%.2f',roc2_1$ci),collapse=('_')))
    result_df$blastpb_ci = c(result_df$blastpb_ci,paste(sprintf('%.2f',roc2_2$ci),collapse=('_')))
    result_df$npm1_ci = c(result_df$npm1_ci,paste(sprintf('%.2f',roc3$ci),collapse=('_')))

    result_df = data.frame(result_df)

    # ==
    text = paste(
        c(paste('ssGSEA score=',sprintf('%.2f',roc1_2$auc),sep=''),
        paste('BLAST.BM=',sprintf('%.2f',roc2_1$auc),sep=''),
        paste('BLAST.PB=',sprintf('%.2f',roc2_2$auc),sep=''),
        paste('NPM1=',sprintf('%.2f',roc3$auc),sep='')), collapse='\n'
                )

    p=ggroc(list(ssGSEA_score_cebpa_genes = roc1_2, blastbm = roc2_1, blastpb = roc2_2, npm1 = roc3))+ 
            annotate('text', x = 0.4, y = 0.10, alpha = .99, colour = "black",label = text,hjust=0,family='serif',fontface =2) +
            theme_bw()+theme(text=element_text(size=15,face='bold',family='serif'))+xlab('Specificity')+ylab('Sensitivity')

    if(plot_flag==TRUE){
        ggsave(sprintf('%s/%s_%s_%s.pdf', dir, filename, auc_ic50, drug), dpi=70, width=7.5, height=5) #
    }
    return(list(
        rocp=p,
        result_df=result_df,
        data=list(res_df_cebpa=res_df_cebpa,
                blast_bm_res_list=blast_bm_res_list,
                blast_pb_res_list=blast_pb_res_list,
                npm1_res_list=npm1_res_list)   
    ))

}

### ssgsea_blast_corr
# exp_mat: gene*sample
# cebpa_set: gene vector
# blastbm_df: sample,Blasts.in.BM/PB
# 
ssgsea_blast_corr <- function(exp_mat,cebpa_set,blastbm_df,y,xlab='ssGSEA score',ylab='BLAST',title='ssgsea_score~blast(Vizome)'){
    df_ssgsea_score = ssgsea_score_gene_set(exp_mat, cebpa_set)

    df_fpr_plot_ssgaea_blast = na.omit(merge(df_ssgsea_score,blastbm_df))
    df_fpr_plot_ssgaea_blast[[y]] = df_fpr_plot_ssgaea_blast[[y]]/100
    p=ggscatter(df_fpr_plot_ssgaea_blast, x='ssGSEA_score', y=y,
                palette = "nature", add = "reg.line", conf.int = TRUE,
                xlab=xlab,ylab=ylab,title=title)+
                stat_cor(label.x = 0.16, label.y=1.05) +
                theme(text=element_text(size=16,face='bold',family='serif'))
    print(p)
    return(p)
}



### ssgsea_NPM1_corr
# exp_mat: gene*sample
# cebpa_set: gene vector
# NPM1_df: sample,NPM1_status
# 
ssgsea_NPM1_corr <- function(exp_mat,cebpa_set,NPM1_df){

    ## ssGSEA打分
    df_ssgsea_score = ssgsea_score_gene_set(exp_mat, cebpa_set)

    df_for_plot_npm1_ssgsea = merge(NPM1_df,df_ssgsea_score)%>%filter(sample%in%ITD_sample_ls)

    p=my_boxplot(df_for_plot_npm1_ssgsea, x = "NPM1_status", y = "ssGSEA_score", 
            label_list = list(labs = 'ssgsea_score~NPM1_status', 
            ylab = "ssgsea_score"), para_list = NA, out_name = NA)


}




########################################################################
########################################################################
if(0){
    ## mut data
    ## mut data
    maf_df = read.csv('../../data/VIZOME/LAML_nature_PDC.maf',sep='\t')
    ITD_sample_ls = as.vector((maf_df %>% filter(Hugo_Symbol=='FLT3'&ITDorNOT=='ITD'))['Tumor_Sample_Barcode'][['Tumor_Sample_Barcode']])
    ITD_sample_ls = str_replace_all(paste('X',ITD_sample_ls,sep=''),'-','.')

    ## exp data
    exp_mat = read_tsv('../../data/VIZOME/nature_aml_log2_fpkm.txt')
    # 处理重复基因，滤出需要的基因
    exp_mat = aggregate(.~ Symbol, exp_mat, mean)
    exp_mat = exp_mat %>% column_to_rownames('Symbol')
    colnames(exp_mat) = str_replace_all(paste('X',colnames(exp_mat),sep=''),'-','.')
    sample_ls = colnames(exp_mat)

    ## drug sen data
    drug_sen_df = read.tsv('../../data/VIZOME/nature_aml_drug_sen.tsv')
    drug_sen_df$lab_id = lapply(drug_sen_df$lab_id, function(x) str_replace_all(paste('X',x,sep=''),'-','.'))
    drug_sen_df$lab_id = as.factor(as.character(drug_sen_df$lab_id))
    drug_ls = drug_sen_df[['inhibitor']][!duplicated(drug_sen_df[['inhibitor']])]


    our_drugs = c('Midostaurin','Quizartinib (AC220)','Crenolanib','Gilteritinib (ASP-2215)')

    ND_gene_dir = './NDgene.txt'
    ND_gene_file = read_tsv(ND_gene_dir)


    # get overlap ND
    MIDO_ls = na.omit(ND_gene_file$MIDO)
    QUIZ_ls = na.omit(ND_gene_file$QUIZ)
    GILTER_ls = na.omit(ND_gene_file$GILTER)
    CRENO_ls = na.omit(ND_gene_file$CRENO)
    ND_gene_list = list(MIDO_ls=MIDO_ls,QUIZ_ls=QUIZ_ls,GILTER_ls=GILTER_ls,CRENO_ls=CRENO_ls)

    ND_overlap_gene = get_plot_overlap(ND_gene_list, 4)



    # get overlap DEG
    get_signif_gene <- function(df, top_num, fc){
        if(fc=='up'){
            sel_gene = ( df[order(df$log2FoldChange,decreasing=TRUE),][1:top_num,]%>%filter(log2FoldChange>0) )$Row.names
        }else if(fc=='down'){
            sel_gene = ( df[order(df$log2FoldChange,decreasing=FALSE),][1:top_num,]%>%filter(log2FoldChange<0) )$Row.names
        }else{
            sel_gene = ( df[order(df$padj,decreasing=FALSE),][1:top_num,] )$Row.names
        }
        return(sel_gene)
    }

    crenolanib_DEG_df = read_csv('/home/lgh/my_project/hl/分训练集测试集/resource/vizome_crenolanib_diff_gene.csv')
    gilterinib_DEG_df = read_csv('/home/lgh/my_project/hl/分训练集测试集/resource/vizome_gilterinib_diff_gene.csv')
    Midostaurin_DEG_df = read_csv('/home/lgh/my_project/hl/分训练集测试集/resource/vizome_Midostaurin_diff_gene.csv')
    quizartinib_DEG_df = read_csv('/home/lgh/my_project/hl/分训练集测试集/resource/vizome_quizartinib_diff_gene.csv')

    crenolanib_signif_gene = get_signif_gene(crenolanib_DEG_df,50,'all')
    gilterinib_signif_gene = get_signif_gene(gilterinib_DEG_df,50,'all')
    Midostaurin_signif_gene = get_signif_gene(Midostaurin_DEG_df,50,'all')
    quizartinib_signif_gene = get_signif_gene(quizartinib_DEG_df,50,'all')

    DEG_gene_list = list(MIDO_ls=Midostaurin_signif_gene,
                        QUIZ_ls=quizartinib_signif_gene,
                        GILTER_ls=gilterinib_signif_gene,
                        CRENO_ls=crenolanib_signif_gene)

    DEG_overlap_gene = get_plot_overlap(DEG_gene_list, 4)





    # test

    drug = our_drugs[1]
    print(drug)
    ## choose sp sample
    itd_exp_drug_sample = intersect(intersect(ITD_sample_ls,colnames(exp_mat)), (drug_sen_df%>%filter(inhibitor==drug))$lab_id)


    # test for ND overlap gene
    test_gene_ls_in_drug(drug_sen_df=drug_sen_df, drug=drug,
                        sp_sample_ls=itd_exp_drug_sample, exp_mat=exp_mat, gene_ls=ND_overlap_gene)

    # test for DEG overlap gene
    test_gene_ls_in_drug(drug_sen_df=drug_sen_df, drug=drug,
                        sp_sample_ls=itd_exp_drug_sample, exp_mat=exp_mat, gene_ls=DEG_overlap_gene)




}
