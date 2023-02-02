# example pipeline: /home/lgh/my_project/omic/modules/grouping/grouping_by_clustering_pipeline.R
# ~~输出的结果需要不同批次结果自动保存成不同的文件夹中
# do_ConsensusCluster
#
#
library(tidyverse)
library(DESeq2)
library(GSVA)
library(pheatmap)
library(ConsensusClusterPlus)
library(ggfortify)
PATH_SRC = '/home/lgh/my_project/omic'
source(sprintf('%s/modules/clean_data/normalization.R',PATH_SRC))


ccResult2annoDf <- function(n_cluster, cc_results, colname='Cluster'){

    class_df = as.data.frame(cc_results[[n_cluster]][['consensusClass']],col.names=c('class'))
    colnames(class_df)=c('class')
    class_df['sample'] = rownames(class_df)

    sp_cc_result = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))

    list_for_factor = c()
    list_for_rn = c()
    for(i in 1:n_cluster){
        list_for_factor = c(list_for_factor, rep(sprintf('c%s',i), length(sp_cc_result[[i]])))
        list_for_rn = c(list_for_rn, as.vector(sp_cc_result[[i]]))
    }

    annotation_df = data.frame(
    colname = factor(list_for_factor)
    )
    rownames(annotation_df) = list_for_rn
    colnames(annotation_df) = colname

    return(annotation_df)

}

make_special_annotation <- function(){
    
}

cc_pheatmap_visualization <- function(exp_mat_cc, 
                                      annotation_col, annotation_row, 
                                      cellwidth, cellheight, 
                                      project_name,result_dir,
                                      colors=c('blue','red'),
                                      other_fn_anno=NULL,fontsize_row=10,fontsize_col=10,
                                      show_rownames=FALSE,show_colnames=FALSE,scale=TRUE){

    # 控制柱饰条颜色
    #anno_colors = list(
    #            Cluster_col = c(c1="red",c2="blue",c3="green",c4="pink",c5='black',c6='yellow',c7='grey'),
    #            Cluster_row = c(c1_g="red",c2_g="blue",c3_g="green",c4_g="pink",c5_g='black'),
    #            term_row = c(CMA='Cyan',UPR='Gold',UPS='Coral',UPS_CMA='DimGray',X='DarkViolet')
    #            )

    # 控制极端值，保证热图好看
    if(scale==TRUE){
        exp_mat_cc_tmp = t(scale(t(exp_mat_cc)))
        bk <- c(seq(-5,-0.1,by=0.01),seq(0,5,by=0.01))
    }else{
        exp_mat_cc_tmp = exp_mat_cc
        bk <- c(seq(min(subtyping_res1$exp_mat),-0.1,by=0.01),seq(0,max(subtyping_res1$exp_mat),by=0.01))
    }

    
    #exp_mat_cc_tmp[exp_mat_cc_tmp>5] = 5
    #exp_mat_cc_tmp[exp_mat_cc_tmp<-5] = -5

    # 图例色条调整


    if(is.null(other_fn_anno)==TRUE){
        filename = paste(c(result_dir,sprintf('cencluster_hmap_ccgene_ccsample_%s.pdf',project_name)), collapse='/')
    }else{
        filename = paste(c(result_dir,sprintf('cencluster_hmap_ccgene_ccsample_%s.pdf',paste(project_name,other_fn_anno,sep='_'))), collapse='/')
    }
    
    pheatmap(
        exp_mat_cc_tmp[rownames(annotation_row),
                   rownames(annotation_col)],
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        #scale = 'row',
        #color = colorRampPalette(c("blue","white", "red"))(50),
        color = c(colorRampPalette(colors = c(colors[1],"white"))(length(bk)/2),colorRampPalette(colors = c("white",colors[2]))(length(bk)/2)),
        #legend_breaks=seq(-5,5,2.5),
        breaks=bk,
        # filename='./cencluster_hmap_ccgene_ccsample_modi_color.pdf',
        filename = filename,
        cellwidth = cellwidth, cellheight = cellheight,
        annotation_col = annotation_col,
        annotation_row = annotation_row,
        #annotation_colors = anno_colors,
        show_rownames=show_rownames,
        show_colnames=show_colnames,
        border_color = NA,
        fontsize_row = fontsize_row,
        fontsize_col = fontsize_col,
        )

}



clean_data_char2num <- function(tmp_exp_mat,log2){
    tmp_exp_mat[tmp_exp_mat=='NA'] = NA
    tmp_exp_mat[,c(2:dim(tmp_exp_mat)[2])] = sapply(tmp_exp_mat[,c(2:dim(tmp_exp_mat)[2])], as.numeric)
    tmp_exp_mat = tmp_exp_mat%>%column_to_rownames('symbols')
    return(tmp_exp_mat)
}

# do cc cluster
#
do_ConsensusCluster <- function(exp_mat_cc,scale=TRUE,
        maxK_sample=8,maxK_gene=8,
        clusterAlg='km',distance='euclidean',
        pItem=0.8,reps=1000, pdf_png='pdf',
        result_dir=NULL,project_name=NULL){
    # exp_mat_cc = gsva_matrix
    exp_mat_cc = na.omit(exp_mat_cc)
    #
    cat(sprintf('## clusterAlg:%s, distance:%s, pItem:%s, reps:%s \n result_dir:%s, project_name:%s\n',
                clusterAlg,distance,pItem,reps,result_dir,project_name))

    # scaling
    if(scale){
        exp_mat_cc = t(scale(t(exp_mat_cc)))
    }
    #

    # 对样本聚类
    # Consensus Cluster #################
    # input of ConsensusClusterPlus: gene×sample ,this func is do the cluster for sample
    ##
    if(!is.null(result_dir)&!is.null(project_name)){

        results_for_sample = ConsensusClusterPlus(exp_mat_cc,maxK=maxK_sample,
                        clusterAlg=clusterAlg,distance=distance,
                        pItem=pItem,reps=reps,
                        title=paste(sprintf('%s/',result_dir),project_name,sep=''),
                        plot=pdf_png) # <- 
    }else{
        results_for_sample = ConsensusClusterPlus(exp_mat_cc,maxK=maxK_sample,
                        clusterAlg=clusterAlg,distance=distance,
                        pItem=pItem,reps=reps) # <-         
    }

    # results = ConsensusClusterPlus(as.matrix(exp_mat_cc),maxK=maxK,clusterAlg='km',distance="Euclidean",pItem=0.8,reps=1000,
    #                              title=paste(sprintf('%s/',work_dir),paste(project_name,pn_top,pna,sep='_'),sep=''),
    #                              plot='pdf')  # 这里设置有讲究



    # 对基因聚类
    # Consensus Cluster #################
    # input of ConsensusClusterPlus: gene×sample ,this func is do the cluster for sample
    ##
    results_for_gene = ConsensusClusterPlus(t(exp_mat_cc),maxK=maxK_gene)

    if(!is.null(result_dir)){
        save(results_for_sample, 
            file = sprintf( paste(c(result_dir,sprintf('consenCluster_res_%s_for_sample.RData',project_name)), collapse='/') ))
        save(results_for_gene, 
            file = sprintf( paste(c(result_dir,sprintf('consenCluster_res_%s_for_gene.RData',project_name)), collapse='/') ))
    }

    subtyping_res = list(
        exp_mat=exp_mat_cc,
        results_for_sample=results_for_sample,
        results_for_gene=results_for_gene
    )

    if(!is.null(result_dir)){        
        save(subtyping_res, 
            file=sprintf( paste(c(result_dir,sprintf('subtyping_res_%s.RData',project_name)), collapse='/') )
            )
    }

    return(subtyping_res)
}

# 获取分类结果df和可视化
## 要解决可视化的label颜色问题
## 解决多层label的问题
#
## colnames_define,rownames_define named vector
#
do_ConsensusCluster_result <- function(subtyping_res,
        n_cluster_sample=3, n_cluster_gene=3,
        rownames_define=NULL,colnames_define=NULL,
        annotation_row_define=NULL,annotation_col_define=NULL,
        result_dir=NULL,project_name=NULL,hmap_label=NULL,
        fontsize_row=10,fontsize_col=NULL,wid_height_scale=1,
        show_rownames=FALSE,show_colnames=FALSE,scale=TRUE){

    post_subtyping_res = list()

    if(is.null(fontsize_col)){
        fontsize_col = dim(subtyping_res$exp_mat)[1]/(8*15/wid_height_scale)
        # print(dim(subtyping_res$exp_mat)[1]/(8*20/wid_height_scale))
    }

    ## make col & row df ##

    # n_cluster_sample = 3
    annotation_col = ccResult2annoDf(n_cluster=n_cluster_sample, 
                                    cc_results=subtyping_res$results_for_sample, 
                                    colname='Cluster_col')

    sample_group_df = data.frame(sample=rownames(annotation_col),group=annotation_col$Cluster_col)
    # save(sample_group_df,file=sprintf('%s/sample_group_df.rda',result_dir))
    # ==============
    # n_cluster_gene = 4
    annotation_row = ccResult2annoDf(n_cluster=n_cluster_gene, 
                                    cc_results=subtyping_res$results_for_gene, 
                                    colname='Cluster_row')
    gene_group_df = data.frame(gene=rownames(annotation_row),group=annotation_row$Cluster_row)
    # save(gene_group_df,file=sprintf('%s/gene_group_df.rda',result_dir))

    post_subtyping_res$sample_group_df = sample_group_df
    post_subtyping_res$gene_group_df = gene_group_df
    post_subtyping_res$result_dir = result_dir
    post_subtyping_res$project_name = project_name

    ## Visualization ##
    # annotation_col
    # annotation_row

    if(is.null(result_dir)|is.null(project_name)){
        return(post_subtyping_res)
    }

    tmp_exp_mat = subtyping_res$exp_mat
    # rownames(tmp_exp_mat) = gene_label
    #pathway_anno_df = read_tsv('./data/pathway_anno_df.txt')
    #pathway_anno_df = pathway_anno_df%>%column_to_rownames('pathway')
    #colnames(pathway_anno_df)[1] = 'Cluster_row'

    #tmp_exp_mat = exp_mat_cc[rownames(pathway_anno_df),]
    if(!is.null(annotation_row_define)){
        annotation_row = annotation_row_define
    }
    if(!is.null(annotation_col_define)){
        annotation_col = annotation_col_define
    }
    tmp_exp_mat = tmp_exp_mat[rownames(annotation_row),rownames(annotation_col)]

    cc_pheatmap_visualization(exp_mat_cc=tmp_exp_mat, 
                            annotation_col=annotation_col, annotation_row=annotation_row, 
                            cellwidth=(dim(tmp_exp_mat)[1]/dim(tmp_exp_mat)[2] * 1.2)*wid_height_scale, cellheight=1*wid_height_scale, 
                            project_name=paste(c(project_name,hmap_label),collapse='_'),
                            result_dir=sprintf('%s',result_dir),fontsize_row=fontsize_row,fontsize_col=fontsize_col,
                            show_rownames=show_rownames,show_colnames=show_colnames,scale=scale)
    marker = 0
    if(!is.null(rownames_define)){
        rownames(tmp_exp_mat) = rownames_define[rownames(tmp_exp_mat)]
        rownames(annotation_row) = colnames_define[rownames(annotation_row)]
        marker = 1
    }
    if(!is.null(colnames_define)){
        colnames(tmp_exp_mat) = colnames_define[colnames(tmp_exp_mat)]
        rownames(annotation_col) = colnames_define[rownames(annotation_col)]
        marker = 1
    }
    if(marker==1){
        cc_pheatmap_visualization(exp_mat_cc=tmp_exp_mat, 
                            annotation_col=annotation_col, annotation_row=annotation_row, 
                            cellwidth=(dim(tmp_exp_mat)[1]/dim(tmp_exp_mat)[2] * 1.2)*20, cellheight=1*20, 
                            project_name=paste(c(project_name,hmap_label,'rowcol'),collapse='_'),
                            result_dir=sprintf('%s',result_dir),fontsize_row=fontsize_row,scale=scale)
    }
    
    write_tsv(sample_group_df,sprintf('%s/%sample_group_df.tsv',result_dir,project_name))
    write_tsv(gene_group_df,sprintf('%s/%s_gene_group_df.tsv',result_dir,project_name))

    save(post_subtyping_res,file=sprintf('%s/post_subtyping_res_%s.rda',result_dir,project_name))

    return(post_subtyping_res)
}



# ======================================
# ======================================

library(FactoMineR)
library(factoextra)
library(rvest)
library(tidyverse)

do_HCPC_cluster <- function(exp_mat_cc,ncp=5,fig=TRUE){
    # Compute PCA with ncp = 3
    res.pca <- PCA(t(exp_mat_cc), ncp = ncp, graph = FALSE)
    # Compute hierarchical clustering on principal components
    res.hcpc <- HCPC(res.pca, graph = FALSE)
    fviz_dend_fig = NULL
    fviz_cluster_fig = NULL
    if(fig==TRUE){
        fviz_dend_fig = fviz_dend(res.hcpc, 
                cex = 0.7,                     # Label size
                palette = "jco",               # Color palette see ?ggpubr::ggpar
                rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
                rect_border = "jco",           # Rectangle color
                labels_track_height = 0.8      # Augment the room for labels
                )

        fviz_cluster_fig = fviz_cluster(res.hcpc,
                    repel = TRUE,            # Avoid label overlapping
                    show.clust.cent = TRUE, # Show cluster centers
                    palette = "jco",         # Color palette see ?ggpubr::ggpar
                    ggtheme = theme_minimal(),
                    main = "Factor map"
                    )
    }
    # plot(res.hcpc, choice = "3D.map")

    return(list(res.pca=res.pca,res.hcpc=res.hcpc,
        fviz_dend_fig=fviz_dend_fig,fviz_cluster_fig=fviz_cluster_fig))

}









#get_rda <- function(rda){
#    l <- load(rda)
#    data <- eval(parse(text = l))
#    return(data)
#}
source('/home/lgh/my_project/omic/modules/utils/load_file.R')

if(FALSE){

    # ==============
    pw_list_path = './data/pw_list.txt'
    exp_path = '~/my_project/data/TCGA/BRCA/TCGA-BRCA.htseq_counts.tsv.gz'
    probeMap_path = '~/my_project/data/TCGA/BRCA/gencode.v22.annotation.gene.probeMap'
    cleaned_exp_path = '...'

    project_name = 'upr_BRCA_only_tumor'
    result_dir = './result/'

    dir.create(result_dir)
    subtyping_res = list()

    ###  create_folder!!!

    ##
    ##

    # === pw list import
    pw_list_df = read_tsv(pw_list_path)

    # ==== TCGA data preprocess

    exp_mat = read_tsv(exp_path)
    probeMap_df = read.csv(probeMap_path,sep='\t')
    ###
    probeMap_df <- probeMap_df %>%  dplyr::select(1,2)
    colnames(probeMap_df) <- c("Ensembl_ID","gene_name")

    exp_mat[,c(2:dim(exp_mat)[2])] = 2**exp_mat[-1]-1 # for count data, reverse

    exp_mat <- inner_join(probeMap_df, exp_mat, by = "Ensembl_ID") %>% dplyr::select(-1)
    exp_mat <- aggregate(.~ gene_name, exp_mat, mean)
    exp_mat <- exp_mat %>% column_to_rownames("gene_name")


    upr_pw_list_df = pw_list_df%>%filter(A=='UPR')



    # ===========
    # ===========

    # =========== TCGA deseq2 normalization

    deseq2_exp_mat = exp_mat
    condition_table = factor(c(rep('sen',round(dim(deseq2_exp_mat)[2]/2)),
                            rep('unsen',dim(deseq2_exp_mat)[2]-round(dim(deseq2_exp_mat)[2]/2))))
    # deseq2_exp_mat = as.matrix(exp_mat%>%column_to_rownames('gene_id'))
    deseq2_exp_mat = apply(deseq2_exp_mat,2,as.integer)
    rownames(deseq2_exp_mat) = rownames(exp_mat)
    deseq2_normalization_res = deseq2_normalization(deseq2_exp_mat,condition_table)

    deseq2_normalized_exp_mat = deseq2_normalization_res$norm_exp_mat

    tumor_sample = as.vector(sapply(substr(colnames(deseq2_normalized_exp_mat),14,15),as.numeric))<10
    deseq2_normalized_exp_mat_bak = deseq2_normalized_exp_mat
    deseq2_normalized_exp_mat = deseq2_normalized_exp_mat[,tumor_sample]




    # ====== clustering ===========
    # pw_genes
    # exp_mat_cc
    pw_genes = unique(upr_pw_list_df$F)

    exp_mat_cc = as.data.frame(deseq2_normalized_exp_mat)[pw_genes,] # 选择pathway基因
    exp_mat_cc = na.omit(exp_mat_cc)
    #
    subtyping_res$exp_mat = deseq2_normalized_exp_mat
    subtyping_res$pw_genes = pw_genes
    # subtyping_res$sel_sample = tumor_sample
    # 对样本聚类
    # Consensus Cluster #################
    # input of ConsensusClusterPlus: gene×sample ,this func is do the cluster for sample
    ##
    results_for_sample = ConsensusClusterPlus(t(scale(t(exp_mat_cc))),maxK=8,reps=1000)
    save(results_for_sample, file = sprintf( paste(c(result_dir,sprintf('consenCluster_res_%s_for_sample.RData',project_name)), collapse='/') ))

    # 对基因聚类
    # Consensus Cluster #################
    # input of ConsensusClusterPlus: gene×sample ,this func is do the cluster for sample
    ##
    results_for_gene = ConsensusClusterPlus(scale(t(exp_mat_cc)),maxK=10,reps=1000)
    save(results_for_gene, file = sprintf( paste(c(result_dir,sprintf('consenCluster_res_%s_for_gene.RData',project_name)), collapse='/') ))

    subtyping_res$results_for_sample = results_for_sample
    subtyping_res$results_for_gene = results_for_gene
    # 


    #######################################
    #######################################
    i_k = 5
    # visualization
    if(FALSE){
        ## xbp1 gene import
        xbp1_genes_table = read_tsv('./data/XBP1_activates_chaperone_genes.txt',col_names=FALSE)
        xbp1_genes = sapply(strsplit(xbp1_genes_table$X3, ' '),'[[',2)
        gene_label = rownames(exp_mat_cc)
        for(i in 1:length(gene_label)){
            if(gene_label[i]%in%xbp1_genes){
                gene_label[i] = paste(gene_label[i],'xbp1')
            }
        }
    }

    n_cluster_sample = i_k
    annotation_col = ccResult2annoDf(n_cluster=n_cluster_sample, 
                                    cc_results=results_for_sample, 
                                    colname='Cluster_col')

    sample_group_df = data.frame(sample=rownames(annotation_col),group=annotation_col$Cluster_col)
    # save(sample_group_df,file=sprintf('%s/sample_group_df.rda',result_dir))
    write_tsv(sample_group_df,file=sprintf('%s/%s_gene_group_df.tsv',result_dir,project_name))

    # ==============
    n_cluster_gene = 5
    annotation_row = ccResult2annoDf(n_cluster=n_cluster_gene, 
                                    cc_results=results_for_gene, 
                                    colname='Cluster_row')
    gene_group_df = data.frame(gene=rownames(annotation_row),group=annotation_row$Cluster_row)
    # save(gene_group_df,file=sprintf('%s/gene_group_df.rda',result_dir))
    write_tsv(gene_group_df,file=sprintf('%s/%s_gene_group_df.tsv',result_dir,project_name))

    subtyping_res$sample_group_df = sample_group_df
    subtyping_res$gene_group_df = gene_group_df
    #
    if(FALSE){
        for(i in 1:length(rownames(annotation_row))){
            if(rownames(annotation_row)[i]%in%xbp1_genes){
                rownames(annotation_row)[i] = paste(rownames(annotation_row)[i],'xbp1')
            }
        }
    }

    tmp_exp_mat = exp_mat_cc
    # rownames(tmp_exp_mat) = gene_label
    cc_pheatmap_visualization(exp_mat_cc=tmp_exp_mat, 
                            annotation_col=annotation_col, annotation_row=annotation_row, 
                            cellwidth=(dim(exp_mat_cc)[1]/dim(exp_mat_cc)[2] * 1.2)*10, cellheight=1*10, 
                            project_name=project_name,
                            result_dir=sprintf('%s/',result_dir))
}

if(FALSE){ # pca

    # =====================
    # 3 pathway PCA & corr
    # =====================

    ## corr between pathway
    prot_list = list()
    for(item in unique(upr_pw_list_df$B)){
        prot_list[[item]] = unique(as.vector((pw_list_df%>%filter(B==item))$F))
    }
    pw_g_list = prot_list
    exp_mat_gsva = deseq2_normalized_exp_mat[unique(pw_list_df$F),] # 选择pathway基因
    exp_mat_gsva = na.omit(exp_mat_gsva)
    gsva_matrix <- gsva(as.matrix(exp_mat_gsva), pw_g_list, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE, parallel.sz=1)



    #library(rgl)
    #library(plotly)
    #attach(mtcars)

    df_for_plot = as.data.frame(t(gsva_matrix))

    # 研究一个函数进行全部比较的方法
    cor.test(df_for_plot$IRE1, df_for_plot$PERK)
    cor.test(df_for_plot$ATF6, df_for_plot$IRE1)
    cor.test(df_for_plot$ATF6, df_for_plot$PERK)

    ## pca
    exp_mat_pca = t(scale(t(exp_mat_cc)))

    df_for_pca = as.data.frame(t(exp_mat_pca[,rownames(annotation_col)]))
    df_for_pca$samples = rownames(df_for_pca)
    df_for_pca$groups = annotation_col$Cluster_col

    pca1<-df_for_pca[,c(1:(dim(df_for_pca)[2]-2))]%>%prcomp()
    autoplot(pca1,data = df_for_pca,col= 'groups',size=2,
            loadings =F,loadings.label = FALSE,
            frame = TRUE,frame.type='norm',
            label = FALSE, label.size = 3
            ) + theme_classic()



    save.image(sprintf('%s/%s.image.rda',result_dir,project_name))
    if(FALSE){
        subtyping_res = list()

        subtyping_res$exp_mat = deseq2_normalized_exp_mat
        subtyping_res$pw_genes = unique(upr_pw_list_df$F)
        subtyping_res$sel_sample = tumor_sample

        subtyping_res$results_for_sample = results_for_sample
        subtyping_res$results_for_gene = results_for_gene


        sample_group_df = data.frame(sample=rownames(annotation_col),
                                    group=annotation_col$Cluster_col)
        gene_group_df = data.frame(gene=rownames(annotation_row),
                                    group=annotation_row$Cluster_row)
        subtyping_res$sample_group_df = sample_group_df
        subtyping_res$gene_group_df = gene_group_df

        save(subtyping_res, file=sprintf('%s/%s.subtyping_res.rda',result_dir,project_name))
    }
    save(subtyping_res, sprintf('%s/%s.subtyping_res.rda',result_dir,project_name))

    ##########
    # 不同数据的分型的比较
    ##########

    brca_gene_group_df = read_tsv('../Proteostasis/result_new_pipeline/gene_group_df.tsv')
    prad_gene_group_df = gene_group_df

    df_for_comp = merge(brca_gene_group_df, prad_gene_group_df)
    colnames(df_for_comp) = c('gene','brca','prad')
    # df_for_comp = df_for_comp%>%filter(brca%in%c('c1','c4')&prad%in%c('c1','c2'))

    library(ggalluvial)
    p = ggplot(data = df_for_comp,
        aes(axis1 = brca, axis2 = prad)) +
    scale_x_discrete(limits = c("brca", "prad"), expand = c(.2, .05)) +
    #xlab("Demographic") +
    geom_alluvium(aes(fill = brca)) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal()
    p

    ggstatsplot::ggbarstats(df_for_comp,  prad, brca,palette = 'Set2',output='plot')

    # save.image(sprintf('%s/%s.image.rda',result_dir,project_name))


}