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
        cellwidth = 30, cellheight = 0.3,
        annotation_col = annotation_col,
        annotation_row = annotation_row,
        #annotation_colors = anno_colors,
        border_color = "grey"
        )

}


result_dir = './result/'
project_name = 'rnaseq_subtyping'

# 数据导入

normalized_exp_mat = read_csv('~/my_project/data/changhai_data/RNAseq整理/AML_rnaseq_normalization_counts.csv')



# 分型

pn_top = 0.15

gene_var = apply(normalized_exp_mat%>%column_to_rownames('gene_symbol'),1,var)
sel_genes = names(gene_var[base::order(gene_var,decreasing=TRUE)][1:round(length(gene_var)*pn_top)])

exp_mat_cc = normalized_exp_mat%>%filter(gene_symbol%in%sel_genes)%>%column_to_rownames('gene_symbol')

####################################################################
# 对样本聚类
# Consensus Cluster #################
# input of ConsensusClusterPlus: gene×sample ,this func is do the cluster for sample
##
results = ConsensusClusterPlus(t(scale(t(exp_mat_cc))),maxK=8,clusterAlg='km',pItem=0.8,reps=1000)
save(results, file = sprintf( paste(c(result_dir,sprintf('consenCluster_res_%s_for_sample.RData',paste(project_name,pn_top,sep='_'))), collapse='/') ))
####

n_cluster_sample = 7
annotation_col = ccResult2annoDf(n_cluster=n_cluster_sample, cc_results=results, colname='Cluster_col')


####################################################################
# 对基因聚类
# Consensus Cluster #################
# input of ConsensusClusterPlus: gene×sample ,this func is do the cluster for sample
##
results_for_gene = ConsensusClusterPlus(scale(t(exp_mat_cc)),maxK=8)
save(results_for_gene, file = sprintf( paste(c(result_dir,sprintf('consenCluster_res_%s_for_gene.RData',paste(project_name,pn_top,sep='_'))), collapse='/') ))
####

n_cluster_gene = 5
annotation_row = ccResult2annoDf(n_cluster=n_cluster_gene, cc_results=results_for_gene, colname='Cluster_row')



# ======= Visualization

other_anno = paste(sprintf('pn%s',pn_top),
                    sprintf('cns%s',n_cluster_sample),
                    sprintf('cng%s',n_cluster_gene),
                    sep='_')

cc_pheatmap_visualization(other_fn_anno=other_anno, exp_mat_cc=exp_mat_cc, 
                          annotation_col=annotation_col, annotation_row=annotation_row, 
                          cellwidth=30, cellheight=0.3, 
                          project_name=project_name,result_dir=result_dir)













