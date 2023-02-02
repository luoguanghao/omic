library(tidyverse)

###########################
###########################
# pair compare
###########################
###########################

library(limma)

# mat_for_limma是mat: gene*sample
# group_df: $sample ; $group 两列
# g2 vs g1 后比前
# FC是两组的mean值相减
# limma需要log的矩阵
bf.do_limma_new <- function(mat_for_limma, group_df, g1, g2, th_sample_pct) {
    cat(sprintf('foldchange %s : %s\n',g2,g1))
    print(table(group_df$group))
    group_df = data.frame(group_df)
    group_df = group_df%>%filter(sample%in%colnames(mat_for_limma))

    group_df$group<-as.character(group_df$group)

    group_df_raw = group_df
    group_df_raw$group = factor(group_df_raw$group, level=c(g1,g2))


    group_df = rbind(group_df%>%filter(group%in%c(g1)),group_df%>%filter(group%in%c(g2)))
    group_df = group_df%>%filter(sample%in%colnames(mat_for_limma))

    mat_for_limma = mat_for_limma[,group_df$sample]

    # == dele many na genes
    g1_sample = (group_df%>%filter(group==g1&sample%in%colnames(mat_for_limma)))$sample
    apply(mat_for_limma[,g1_sample],1,function(x) sum(!is.na(x)))/length(g1_sample)

    g2_sample = (group_df%>%filter(group==g2&sample%in%colnames(mat_for_limma)))$sample
    apply(mat_for_limma[,g2_sample],1,function(x) sum(!is.na(x)))/length(g2_sample)

    pct_g1 = (apply(mat_for_limma[,g1_sample],1,function(x) sum(!is.na(x)))/length(g1_sample))
    pct_g2 = (apply(mat_for_limma[,g2_sample],1,function(x) sum(!is.na(x)))/length(g2_sample))
    
    mat_for_limma =  mat_for_limma[(pct_g1>=th_sample_pct & pct_g2>=th_sample_pct),]
    # ==

    group_df[group_df$group==g2,]$group = 2
    group_df[group_df$group==g1,]$group = 1
    condition_table = group_df$group
    
    
    design <- model.matrix(~condition_table)
    colnames(design) <- levels(condition_table)
    rownames(design) <- colnames(mat_for_limma)
    fit <- lmFit(mat_for_limma, design)
    fit <- eBayes(fit, trend=TRUE)
    limma_res <- topTable(fit, coef=2,n=Inf)

    limma_res = cbind(limma_res,mat_for_limma[rownames(limma_res),])
    limma_res=limma_res%>%rownames_to_column('row.names')
    colnames(limma_res)[1:7] = c('row.names','log2FoldChange','AveExpr','stat','p','padj','B')
    
    print('## finish limma!')
    return(list(resdata=limma_res, group_df=group_df_raw))

}

###########################
###########################
# multi level compare
###########################
###########################
# compare is tttmp
bf.do_limma_multi_level <- function(mat_for_limma, group_df, compare, th_sample_pct){


    group_df = group_df%>%filter(sample%in%colnames(mat_for_limma))
    mat_for_limma = (mat_for_limma[,group_df$sample])

    # == dele many na genes
        
    group_vec = unique(group_df$group)
    for(g in group_vec){
        g1_sample = (group_df%>%filter(group==g&sample%in%colnames(mat_for_limma)))$sample
        apply(mat_for_limma[,g1_sample],1,function(x) sum(!is.na(x)))/length(g1_sample)
              
        pct_g1 = (apply(mat_for_limma[,g1_sample],1,function(x) sum(!is.na(x)))/length(g1_sample))
        
        mat_for_limma =  mat_for_limma[(pct_g1>=th_sample_pct),]
    }                
    
    # ==
    
    
    group = group_df$group
    design <- model.matrix(~0+factor(group))
    colnames(design)=levels(factor(group))
    rownames(design)=colnames(mat_for_limma)

    contrast.matrix = eval(parse(text = 
        sprintf('makeContrasts(%s,levels=design)',paste(sprintf('"%s"',compare),collapse=','))))

    ##step1
    fit <- lmFit(mat_for_limma,design)
    ##step2
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)  
    ##step3
    tempOutput = topTable(fit2, coef=1, n=Inf)
    nrDEG = na.omit(tempOutput) 

    result = list(overall_test_result = topTable(fit2, coef=NULL, n=Inf)%>%rownames_to_column('gene'))
    for(i_comp in 1:length(compare)){
        result[[ compare[i_comp] ]] = topTable(fit2, coef=i_comp, n=Inf)%>%rownames_to_column('gene')
        colnames(result[[ compare[i_comp] ]]) = c('gene','log2FoldChange','AveExpr','t','p','padj','B')
        result[[ compare[i_comp] ]] = list(resdata=result[[ compare[i_comp] ]])
    }

    # ==add mean==
    exp_df = mat_for_limma%>%rownames_to_column('gene')
    colnames(exp_df)[1] = 'gene'
    data3 <- reshape2::melt(exp_df, id.vars = c("gene"))
    colnames(data3) = c('variable','sample','value')

    #group_df=as.data.frame(immu_subtype_df%>%filter(group%in%c('C1','C2','C3')))

    df_for_cal_pro = merge(group_df,data3)%>%filter(!is.na(value))

    group_mean_pro = df_for_cal_pro%>%group_by(variable,group)%>%summarize(mean=mean(value))
    tmp_mean = reshape2::dcast(group_mean_pro, variable~group, value.var='mean')
    colnames(tmp_mean)[2:4] = paste('mean_',colnames(tmp_mean)[2:4],sep='')
    colnames(tmp_mean)[1] = 'gene'
                        
    result$overall_test_result = merge(result$overall_test_result,tmp_mean,by='gene')
    # ============

    return(result)

}

bf.do_limma_multi_level.make_comb_df <- function(do_limma_multi_level_res, adj=FALSE){
    result_df = do_limma_multi_level_res[['overall_test_result']]
    for(i_comp in 2:length(names(do_limma_multi_level_res))){
        if(adj==FALSE){
            tmp_df = do_limma_multi_level_res[[i_comp]][,c('gene','logFC','P.Value')]

            colnames(tmp_df)[2:3] = paste(colnames(result_df)[i_comp],colnames(tmp_df)[2:3],sep='.')
            result_df = merge(result_df,tmp_df,by='gene')
        }else{ # 
            tmp_df = do_limma_multi_level_res[[i_comp]][,c('gene','logFC','P.Value','adj.P.Val')]

            colnames(tmp_df)[2:4] = paste(colnames(result_df)[i_comp],colnames(tmp_df)[2:4],sep='.')
            result_df = merge(result_df,tmp_df,by='gene')
        }

    }
    result_df = result_df[,-c(2:length(names(do_limma_multi_level_res)))]
    return(result_df)
}


###########################
###########################
# multi level compare
###########################
###########################
# compare is tttmp
# normalized_exp_mat: gene*sample

bf.prepare_gsea <- function(normalized_exp_mat, group_df, g1, g2, gmtfile, result_dir, comp_name, num_show){
  print(sprintf('## %s/%s_exp_mat.gct',result_dir,comp_name))
  
  dir.create(result_dir,recursive = T)

  out = result_dir

  rownames(group_df) = NULL

  common_samples = intersect(group_df$sample, colnames(normalized_exp_mat))
  condition_table = (data.frame(group_df)%>%column_to_rownames('sample'))[common_samples,]
  normalized_exp_mat = normalized_exp_mat[,common_samples]
  
  prepare_gsea_file(normalized_exp_mat=normalized_exp_mat,
        condition_table=condition_table,
        result_dir=result_dir,
        comp_name=comp_name)

  cmd_list = list()

  for(i in 1:length(gmtfile)){
    gmt_name = names(gmtfile)[i]
    gmt=gmtfile[[i]]

    gct=sprintf('"%s/%s_exp_mat.gct"',result_dir,comp_name)
    cls=sprintf('"%s/%s.cls"',result_dir,comp_name)
    #g1='M12'
    #g2='M345'
    rpt_label=sprintf('"%s_%s"',comp_name,gmt_name)
    metric='Signal2Noise'
    permute='geneset'
    # out='./fab/'
    cmd = gsea_cmd( gct, cls, g1, g2, gmt, rpt_label, metric, permute, sprintf('"%s"',out), num_show=num_show )
    #print(gmt_name)
    #print(cmd)
    cmd_list[[gmt_name]] = cmd
  }


  return(cmd_list)

}





