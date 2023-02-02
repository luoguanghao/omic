source('/home/lgh/my_project/omic/modules/feature_selection/cor_one_by_one.R')
source('/home/lgh/my_project/omic/modules/bioinformatics_wrapper/enrich.R')
source('/home/lgh/my_project/omic/modules/grouping/func_post_subtyping.R')
library(tidyverse)
options ( warn = -1)

############################
# wilcox
############################

compare_proteome <- function(protein_df, group_df, gmtfile_ls,  method='wilcox',g1,g2, log=FALSE,gsea=FALSE){
    
    diff_exp_one_by_one_res = diff_exp_one_by_one(exp_df=protein_df, group_df=group_df, 
                                                  g1=g1,g2=g2, method='wilcox', 
                                                  genes=NULL, t=FALSE, log=log)
    
    # enrichment
    gsea_res = list()
    if(gsea==TRUE){ 
        for(gsea_type in names(gmtfile_ls)){
        gene_df = data.frame(gene=diff_exp_one_by_one_res$test_df$gene, score_for_rank=diff_exp_one_by_one_res$test_df$log2_fc)  ###
        #gmtfile = paste(PATH_DATA,'/pathway/c5.go.bp.v7.2.symbols.gmt',sep=''  )            ### go all
        gsea_res[[gsea_type]] = GSEA_enrich(gene_df, gmtfile_ls[[gsea_type]], pvalueCutoff=1)      
        }
    }

    return(list(test_df=diff_exp_one_by_one_res$test_df,
                diff_exp_one_by_one_res=diff_exp_one_by_one_res,
                gsea_res=gsea_res
               ))
}

compare_phoshpoproteome <- function(phosph_df, KSEAdata=NULL, group_df,  method='wilcox',g1,g2, log=FALSE){
    
    diff_exp_one_by_one_res = diff_exp_one_by_one(exp_df=phosph_df, group_df=group_df, 
                                                g1=g1,g2=g2, method='wilcox', 
                                                genes=NULL, t=FALSE, log=log)
    # return(diff_exp_one_by_one_res)
    ####
    # anno kinase
    anno_kinase_res = anno_kinase(diff_exp_one_by_one_res$test_df,diff_exp_one_by_one_res$test_df_with_count)
    diff_exp_one_by_one_res$test_df_anno_kinase = anno_kinase_res$test_df_anno_kinase
    diff_exp_one_by_one_res$test_df_with_count_anno_kinase = anno_kinase_res$test_df_with_count_anno_kinase
    ####
    ####

    # KSEA
    ## effect_name,p_name='p',t=FALSE,unlog=FALSE
    #do_KSEA_res = do_KSEA(pho_df=phosph_df,
    #        stat_res=diff_exp_one_by_one_res$test_df,
    #        effect_name='log2_fc',t=TRUE,unlog=FALSE)    
    #ksea_scores = do_KSEA_res$ksea_scores[order(do_KSEA_res$ksea_scores$p.value),]
    
    return(list(test_df=diff_exp_one_by_one_res$test_df,
               diff_exp_one_by_one_res=diff_exp_one_by_one_res,
               ksea_scores=NA,
               do_KSEA_res=NA))
}

############################
# survival
############################







############################
# plot
############################
# stat_df: point_name, r, p, p.adj, other...
# sp_points = c('BCL2','BCL2A1','MCL1')
# effect_name = 'r' or 'fc'
# label_name = 'gene' ...
plot_volcano_old <- function(stat_df,sp_points,effect_name,label_name,title,thr=0.4){

    stat_df$effect = stat_df[[effect_name]]
    stat_df$label = stat_df[[label_name]]

    volc = ggplot(stat_df, aes(x=`effect`, y=-log10(`p.adj`))) +
        geom_point() + # scale_color_manual(values=c("#4169E1","grey", "#DC143C")) + # 这个和factor顺序有关
        ggtitle(title)+# geom_point(data=data[data$gene %in% marked_genes,], aes(log2FoldChange, -log10(padj)), colour="blue", size=2) +
        geom_hline(yintercept = -log10(0.1),lty=4,lwd=0.6,alpha=0.8)+
        geom_vline(xintercept = c(0.4,-0.4),lty=4,lwd=0.6,alpha=0.8)

    #volc = volc+ggrepel::geom_text_repel(data=drugsen_gene_res[drugsen_gene_res[['gene']]%in%c('BCL2','BCL2A1','MCL1'),], aes(label=`gene`), max.overlaps=20) + theme_bw(base_size=15)
    # 可以换成 volc = volc+ggrepel::geom_label_repel(data=df_for_plot[df_for_plot[['Cancer feature']]=='FLT3_mut',], aes(label=`Cancer feature`), max.overlaps=20,size=5)
    volc = volc + geom_point(data=stat_df[stat_df[[label_name]]%in%sp_points,], 
                            aes(x=`effect`, y=-log10(`p.adj`)), color='red',size=5)
    volc = volc+ggrepel::geom_text_repel(data=stat_df[stat_df[[label_name]]%in%sp_points,], 
                            aes(label=`label`), max.overlaps=20, color='red') + theme_bw(base_size=15)

    plot(volc)

    sp_df = stat_df%>%filter(label%in%sp_points)

    return(list(plot=volc,sp_df=sp_df))
}



###############################
# wrapped diff & diff & fisher
###############################
# 样本矩阵注意列名以字母开头！！
#
# gsea_data_list=list(normalized_exp_mat=normalized_exp_mat, 
#                     result_dir='./',comp_name='Pair_GLS_vs_EtOH_GLS_0507',gmtfile=gmtfile, num_show =200)
#
omics_diff_analysis <- function(group_df=NULL,g1,g2,
                        suv_data_list=NULL,exp_data_list=NULL,
                        protein_data_list=NULL, phosph_data_list=NULL,phosphgene_data_list=NULL,
                        mut_data_list=NULL,drug_data_list=NULL,
                        gsea_data_list=NULL,aucic50='ic50'){
    
    post_subtyping_res = list()
    
    # survival
    if(!is.null(suv_data_list)){

        # compare_survival_res = compare_survival(suv_data_list$suv_df, group_df=group_df%>%filter(group%in%c('c1','c2','c3')),plot=FALSE)
        compare_survival_res = compare_survival(suv_data_list$suv_df, group_df=group_df,plot=suv_data_list$plot)
        post_subtyping_res$suv = compare_survival_res  ## <--

        #if(suv_data_list$pairwise==TRUE){   
        #}

    }    

    # exp
    if(!is.null(exp_data_list)){ # add multi_compare and pairwise_compare
        
        compare_transcriptome_res = compare_transcriptome_pair(exp_mat=exp_data_list$exp_mat,
                gmtfile=exp_data_list$gmtfile,
                group_df=group_df, g1=g1, g2=g2, ## <--
                method=exp_data_list$method,gsea = FALSE)
        post_subtyping_res$exp = compare_transcriptome_res

        ##
        # make GSEA file
        # prepare_gsea_file(normalized_exp_mat,condition_table,result_dir,comp_name)
        
        ##
        
    }   
   
    # protein
    if(!is.null(protein_data_list)){ # add multi_compare and pairwise_compare
        
        compare_proteome_res = compare_proteome(protein_data_list$protein_df, group_df, 
                exp_data_list$gmtfile,  method='wilcox', g1=g1, g2=g2, log=protein_data_list$log)
        post_subtyping_res$protein = compare_proteome_res

    }      
   
    # phosphoprotein
    if(!is.null(phosph_data_list)){ # add multi_compare and pairwise_compare
        
        compare_phoshpo_res = compare_phoshpoproteome(phosph_data_list$phosph_df, KSEAdata=NULL, group_df,  
                method='wilcox', g1=g1, g2=g2, log=phosph_data_list$log)
        post_subtyping_res$phosph = compare_phoshpo_res

    }

    if(!is.null(phosphgene_data_list)){ # add multi_compare and pairwise_compare
        
        compare_phosphgene_res = compare_proteome(phosphgene_data_list$phosphgene_df, group_df, 
                exp_data_list$gmtfile,  method='wilcox', g1=g1, g2=g2, log=phosphgene_data_list$log)
        post_subtyping_res$phosphgene = compare_phosphgene_res

    }    

    # compare_drug_sen
    # 调 group_df 的 factor
    if(!is.null(drug_data_list)){
        compare_drug_sen_res = compare_drug_sen(drug_data_list$drug_sen_df, group_df, 
                    drug_anno_df=NULL,
                    ic50auc=drug_data_list$ic50auc)

        post_subtyping_res$drug = compare_drug_sen_res
    }

    # prepare_GSEA
    if(!is.null(gsea_data_list)){
        # post_subtyping.prepare_gsea(normalized_exp_mat, group_df, g1, g2, out, gmtfile, result_dir, comp_name)
        prepare_gsea_res = post_subtyping.prepare_gsea(normalized_exp_mat=gsea_data_list$normalized_exp_mat,
                group_df=group_df,
                g1=g1, g2=g2,
                out=gsea_data_list$result_dir,
                gmtfile=gsea_data_list$gmtfile,
                result_dir=gsea_data_list$result_dir,
                comp_name=gsea_data_list$comp_name, num_show=gsea_data_list$num_show)
        
        post_subtyping_res$gsea = prepare_gsea_res
    }



    return(post_subtyping_res)
}


if(FALSE){


    mut_mat_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_mutations_matrix_cleaned.tsv')
    mut_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_mutations_cleaned.maf')

    exp_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_count_expression_cleaned.tsv')
    norm_exp_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_fpkm_expression_cleaned.tsv')

    protein_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_protein_withna_log.tsv')
    phosph_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_phosphoprotein_withna_log.tsv')


    suv_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_survival_fab_cleaned.tsv')

    drug_sen_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_rawic50_drug_sensitivity_cleaned.tsv')

    gmtfile = list(
        go=paste(PATH_DATA,'/pathway/c5.go.bp.v7.2.symbols.gmt',sep=''),
        kegg=paste(PATH_DATA,'/pathway/c2.cp.kegg.v7.2.symbols.gmt',sep=''),
        reactome=paste(PATH_DATA,'/pathway/reactome.gmt',sep='')
    )

    # =====
    suv_data_list=list(suv_df=suv_df,plot=TRUE)
    exp_data_list=list(exp_mat=exp_df%>%column_to_rownames('Hugo_Symbol'), gmtfile=gmtfile, method='deseq2')
    protein_data_list=list(protein_df=protein_df)
    phosph_data_list=list(phosph_df=phosph_df)
    drug_data_list=list(drug_sen_df=drug_sen_df, ic50auc='ic50c')
    gsea_data_list=list(normalized_exp_mat=normalized_exp_mat, 
                        condition_table=condition_table,
                        result_dir,comp_name)

    # =====

    group_df = mut_mat_df[,c('sample','FLT3')]
    colnames(group_df)[2] = 'group'
    group_df$group=as.character(group_df$group)
    group_df[group_df$group==1,]$group='mut'
    group_df[group_df$group==0,]$group='wt'

    omics_diff_analysis_res = omics_diff_analysis(group_df=group_df,
                        suv_data_list=NULL,exp_data_list=NULL,
                        protein_data_list=NULL, phosph_data_list=phosph_data_list,
                        mut_data_list=NULL,drug_data_list=drug_data_list,
                        gsea_data_list=NULL)


}

# exp_df: gene*sample
# t_protein_df,t_pho_df: sample*gene: phosph_site: GENE_AXXX
# drug: just drug name
omics_cor_analysis <- function(drug_sen_df, exp_df=NULL, t_protein_df=NULL, 
                            t_pho_gene_df=NULL, t_pho_df=NULL, drug, aucic50='ic50',method='spearman',orderby='p'){

    all_res_list = list()

    # exp
    if(!is.null(exp_df)){
        drugsen_gene_res = drugsen_gene_cor(drug_sen_df,
            exp_df,drug=drug,
            genes=NULL, t=TRUE, aucic50=aucic50,method=method,orderby=orderby)
            #genes=get_most_var_feature(exp_df%>%column_to_rownames('Hugo_Symbol'),6000))
        
        all_res_list$exp = drugsen_gene_res
    }

    # pro
    if(!is.null(t_protein_df)){
        drugsen_gene_cor_res = drugsen_gene_cor(drug_sen_df=drug_sen_df,
            exp_df=t_protein_df,
            drug=drug,t = FALSE, aucic50=aucic50,method=method,orderby=orderby)

        all_res_list$pro = drugsen_gene_cor_res
    }

    # phosph_gene
    if(!is.null(t_pho_gene_df)){
        drugsen_gene_cor_res = drugsen_gene_cor(drug_sen_df=drug_sen_df,
            exp_df=t_pho_gene_df,
            drug=drug,t = FALSE, aucic50=aucic50,method=method,orderby=orderby)

        all_res_list$phosph_gene = drugsen_gene_cor_res

    }

    # phosph
    if(!is.null(t_pho_df)){
        drugsen_gene_cor_res = drugsen_gene_cor(drug_sen_df=drug_sen_df,
            exp_df=t_pho_df,
            drug=drug,t = FALSE, aucic50=aucic50,method=method,orderby=orderby)

        all_res_list$phosph = drugsen_gene_cor_res

        # KSEA
        do_KSEA_res = do_KSEA(t_pho_df,
                        drugsen_gene_cor_res$drugsen_gene_corr_df,
                        effect_name='r',t=TRUE,unlog = TRUE)
        do_KSEA_res$ksea_scores = do_KSEA_res$ksea_scores[order(do_KSEA_res$ksea_scores$p.value),]
        all_res_list$ksea = do_KSEA_res

    }


    # GSEA


    return(all_res_list)

}



###############################
# 结果输出
###############################


###############################
# 结果分析
###############################

# drug corr 

# plot #
# stat_df: point_name, r, p, p.adj, other...
# sp_points = c('BCL2','BCL2A1','MCL1')
# effect_name = 'r' or 'fc'
# label_name = 'gene' ...
plot_volcano <- function(df_for_plot,drug,omic,effect_size,p,show_lable,show_vec,th_p,th_fc){
    options(repr.plot.width=8, repr.plot.height=8)

    # df_for_plot = omics_cor_analysis_res[['Arac']][['exp']][['drugsen_gene_corr_df']]
    
    df_for_plot$effect_size = df_for_plot[[effect_size]]
    df_for_plot$p = df_for_plot[[p]]
    df_for_plot$show_lable = df_for_plot[[show_lable]]

    volc = ggplot(df_for_plot, aes(x=`effect_size`, y=-log10(p))) +
        geom_point() + scale_color_manual(values=c("#4169E1","grey", "#DC143C")) + # 这个和factor顺序有关
        ggtitle(sprintf("%s %s",drug,omic))+# geom_point(data=data[data$gene %in% marked_genes,], aes(log2FoldChange, -log10(padj)), colour="blue", size=2) +
        geom_hline(yintercept = -log10(th_p),lty=4,lwd=0.6,alpha=0.8)+
        geom_vline(xintercept = c(th_fc,-th_fc),lty=4,lwd=0.6,alpha=0.8)+
        xlab(effect_size)+ylab(sprintf('log10(%s)',p))

    # 可以换成 volc = volc+ggrepel::geom_label_repel(data=df_for_plot[df_for_plot[['Cancer feature']]=='FLT3_mut',], aes(label=`Cancer feature`), max.overlaps=20,size=5)
    volc = volc + geom_point(data=df_for_plot[df_for_plot[['show_lable']]%in%show_vec,], 
                    aes(x=`effect_size`, y=-log10(`p`)), color='red',size=5)
    volc = volc+geom_text_repel(data=df_for_plot[df_for_plot[['show_lable']]%in%show_vec,], 
                    aes(label=show_lable), max.overlaps=20) + theme_bw(base_size=15)

    # ggsave('Wee1_plot.pdf', dpi=50, width=8, height=8)

    return(volc)

}

description_stat <- function(df,effect_size,p,th_fc,th_p){
    df$effect_size = df[[effect_size]]
    df$p = df[[p]]
    cat(
        sprintf('卡 |%s|>%s, %s<%s 的情况下，正有%s个基因，负有%s个基因\n',
               effect_size,th_fc,p,th_p,
               dim(df%>%filter(effect_size>th_fc&p<th_p))[1],
               dim(df%>%filter(effect_size < -th_fc&p<th_p))[1])
    )
}

#
#
kinase_annotation <- function(drugsen_gene_corr_df, kinase_anno, drug, result_dir, effect_size, stat_signif){

    drugsen_gene_corr_df[['effect_size']] = drugsen_gene_corr_df[[effect_size]]
    drugsen_gene_corr_df[['stat_signif']] = drugsen_gene_corr_df[[stat_signif]]

    drugsen_gene_corr_df_kinase_anno = merge(kinase_anno,drugsen_gene_corr_df,
        by.x='gene',
        by.y='gene',all.y=TRUE)

    tmp_df = drugsen_gene_corr_df_kinase_anno[
        order(drugsen_gene_corr_df_kinase_anno$effect_size,decreasing=TRUE)
        ,]%>%filter(stat_signif<0.05&!is.na(kinase))

    tmp_df_mean = tmp_df%>%group_by(kinase)%>%summarise(mean = mean(effect_size))
    tmp_df_mean = tmp_df_mean[order(tmp_df_mean$mean,decreasing=TRUE),]

    #tmp_df_mean$kinase
    tmp_df$kinase = factor(tmp_df$kinase, level=tmp_df_mean$kinase)

    # ===

    # options(repr.plot.width=20, repr.plot.height=8)
    plot = ggplot(tmp_df,aes(x=kinase,y=effect_size)) + geom_point()
    # p = p+geom_text_repel(data=tmp_df, aes(label=`gene`), max.overlaps=20)+theme_bw(base_size=15)
    plot = plot+theme(text=element_text(size=15, face = "bold"),
                    axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.3),
                    panel.background=element_blank(),
                    axis.line = element_line(colour="black"))+labs(title=drug)+
                    geom_hline(yintercept = 0, alpha=0.5, linetype="dashed")
    ggsave(sprintf('%s/%s_phosph_kinase.pdf', result_dir, drug),width=15,height=6)
    # print(p)
    # options(repr.plot.width=8, repr.plot.height=8)

    # ===

    write_tsv(drugsen_gene_corr_df_kinase_anno,
            file=sprintf('%s/%s_drugsen_gene_corr_df_kinase_anno.tsv',result_dir,drug))


    return(list(
        drugsen_gene_corr_df_kinase_anno=drugsen_gene_corr_df_kinase_anno,
        plot=plot     
                )
            )

}


plot_kinase_and_substrate <- function(kinase_anno_df,kinase_ls){
    
    tmp_df = kinase_anno_df[
        order(kinase_anno_df$effect_size,decreasing=TRUE),]%>%filter(kinase%in%kinase_ls)

    tmp_df_mean = tmp_df%>%group_by(kinase)%>%summarise(mean = mean(effect_size))
    tmp_df_mean = tmp_df_mean[order(tmp_df_mean$mean,decreasing=TRUE),]

    #tmp_df_mean$kinase
    tmp_df$kinase = factor(tmp_df$kinase, level=tmp_df_mean$kinase)

    # ===

    # options(repr.plot.width=20, repr.plot.height=8)
    plot = ggplot(tmp_df,aes(x=kinase,y=effect_size)) + geom_point(aes(color=p<0.05),shape='_',size=8)
    # p = p+geom_text_repel(data=tmp_df, aes(label=`gene`), max.overlaps=20)+theme_bw(base_size=15)
    plot = plot+theme(text=element_text(size=15, face = "bold"),
                    axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.3),
                    panel.background=element_blank(),
                    axis.line = element_line(colour="black"))+labs(title=drug)+
                    geom_hline(yintercept = 0, alpha=0.5, linetype="dashed")
    # ggsave(sprintf('%s/%s_phosph_kinase.pdf', result_dir, drug),width=15,height=6)    
    return(list(plot=plot,df_for_plot=tmp_df))
}


compare_drug_signif_genes <- function(res_list,drug1,drug2,p1,p2,omic,df_name,stat_signif,show_label){
    
    df1 = res_list[[drug1]][[omic]][[df_name]]
    df2 = res_list[[drug2]][[omic]][[df_name]]
    df1[['stat_signif']] = df1[[stat_signif]]
    df2[['stat_signif']] = df2[[stat_signif]]

    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    x = list()
    x[[drug1]] = (df1%>%filter(stat_signif<p1))[[show_label]]
    x[[drug2]] = (df2%>%filter(stat_signif<p2))[[show_label]]
    venn.plot = venn.diagram(x,
                filename=NULL,main=sprintf('%s.vs.%s %s',drug1,drug2,omic),main.cex = 2,
                cex = 2.5,cat.cex = 2.5 #,height = 1500, width = 3000
            )
    grid.draw(venn.plot);
    grid.newpage();

    return(get.venn.partitions(x))
}

















