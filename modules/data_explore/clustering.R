library(stringr)
library(ggpubr) 
library(ggplot2)
library(GSEABase)
# library(GSVAdata)
#data(c2BroadSets)
#c2BroadSets

library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)
library(tidyverse)

############
# hcluster #
############
library(pheatmap)
#
# input: df_for_plot, anno_nc_df
## df_for_plot:matrix, context is numeric;colname&rownames
#
# return: df,describe the group of each item in rownames.
# ! only do clustering for columns!!
dohcluster <- function(df_for_plot, scale = 'column', cutree_cols=4, output_file=NULL, param_ls=NULL, test_ls=NULL) {

  p_res = pheatmap(
      df_for_plot,
      #cluster_cols = FALSE,
      cluster_rows = FALSE,
      scale = 'column',
      color = colorRampPalette(c("green","black", "red"))(50),
      # cellwidth = 8, cellheight = 28,border=FALSE,
      cutree_cols=cutree_cols, 
      # annotation_col  = anno_nc_df,
      # main = drug
      # filename='./hmap2.pdf'
      )

  col_cluster <- cutree(p_res$tree_col,k=cutree_cols)

  # get item names
  ## labels(col_cluster[col_cluster==3])

  return(as.data.frame(col_cluster))

}


library(tidyverse)
library(ConsensusClusterPlus)
library(pheatmap)



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


cc_pheatmap_visualization <- function(other_fn_anno, exp_mat_cc, 
                                      annotation_col, annotation_row, 
                                      cellwidth, cellheight, 
                                      project_name,result_dir){

    # 控制柱饰条颜色
    #anno_colors = list(
    #            Cluster_col = c(c1="red",c2="blue",c3="green",c4="pink",c5='black',c6='yellow',c7='grey'),
    #            Cluster_row = c(c1_g="red",c2_g="blue",c3_g="green",c4_g="pink",c5_g='black'),
    #            term_row = c(CMA='Cyan',UPR='Gold',UPS='Coral',UPS_CMA='DimGray',X='DarkViolet')
    #            )

    # 控制极端值，保证热图好看
    exp_mat_cc_tmp = t(scale(t(exp_mat_cc)))
    #exp_mat_cc_tmp[exp_mat_cc_tmp>5] = 5
    #exp_mat_cc_tmp[exp_mat_cc_tmp<-5] = -5

    # 图例色条调整
    bk <- c(seq(-5,-0.1,by=0.01),seq(0,5,by=0.01))


    pheatmap(
        exp_mat_cc_tmp[rownames(annotation_row),
                   rownames(annotation_col)],
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        #scale = 'row',
        #color = colorRampPalette(c("blue","white", "red"))(50),
        color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
        #legend_breaks=seq(-5,5,2.5),
        breaks=bk,
        # filename='./cencluster_hmap_ccgene_ccsample_modi_color.pdf',
        filename = paste(c(result_dir,sprintf('cencluster_hmap_ccgene_ccsample_%s.pdf',paste(project_name,other_anno,sep='_'))), collapse='/'),
        cellwidth = cellwidth, cellheight = cellheight,
        annotation_col = annotation_col,
        annotation_row = annotation_row,
        #annotation_colors = anno_colors,
        border_color = "grey"
        )

}


if(FALSE){
  ####################
  # ConsensusCluster #
  ####################
  # input matrix_cc: feature × sample with rnames & colnames
  library(ConsensusClusterPlus)
  library(pheatmap)
  ###



  # results = ConsensusClusterPlus(t(scale(t(matrix_cc))),maxK=8)
  results = ConsensusClusterPlus(mat), maxK=8)

  ### 这里是取4类的结果来分析,观察result里面给出的最佳类数
  class_df = as.data.frame(results[[4]][['consensusClass']],col.names=c('class'))
  colnames(class_df)=c('class')
  class_df['sample'] = rownames(class_df)

  tmp1 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['1']]
  tmp2 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['2']]
  tmp3 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['3']]
  tmp4 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['4']]

  ### 画热图
  #### 4c
  annotation_col = data.frame(
    Cluster = factor(c(rep('c1',length(tmp1)) , rep('c2',length(tmp2)) , rep('c3',length(tmp3)) , rep('c4',length(tmp4))))
  )
  rownames(annotation_col) = c(tmp1,tmp2,tmp3,tmp4)

  rownames(metabol_df_cc) = toupper(as.character(rownames(metabol_df_cc)))
  pheatmap(
      metabol_df_cc[,c(tmp1,tmp2,tmp3,tmp4)],
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      scale = 'row',
      color = colorRampPalette(c("green","black", "red"))(50),
      # filename='./sorted_cencluster_hmap2.pdf',
      cellwidth = 17, cellheight = 8,
      annotation_col = annotation_col,
      border=FALSE
      )




}





