source(sprintf('%s/modules/clean_data/normalization.R',PATH_SRC))

library(tidyverse)
library(DESeq2)
library(GSVA)
library(pheatmap)
library(ConsensusClusterPlus)


deseq2_normalization <- function(deseq2_exp_mat, condition_table){

    # =========
    dds <- DESeqDataSetFromMatrix(deseq2_exp_mat, 
                                DataFrame(condition_table), 
                                design= ~ condition_table)
    dds <- dds[rowSums(counts(dds)) > 1,]
    dds2 <- DESeq(dds)
    resultsNames(dds2)
    # acquire the results using function results(), and assign to res
    res <- results(dds2)
    # view the summary of results
    summary(res)
    # resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
    # output
    #head(resdata)

    deseq2_normalized_exp_mat = as.data.frame(counts(dds2,normalize=TRUE))

    return(deseq2_normalized_exp_mat)
}


pw_list_path = './data/pw_list.txt'
exp_path = '~/my_project/data/TCGA/BRCA/TCGA-BRCA.htseq_counts.tsv.gz'
probeMap_path = '~/my_project/data/TCGA/BRCA/gencode.v22.annotation.gene.probeMap'
cleaned_exp_path = '...'
result_dir = './result/'
project_name = 'BCRA'

# === pw list import
pw_list_df = read_tsv(pw_list_path)

# ==== TCGA data preprocess

exp_mat = read_tsv(exp_path)
probeMap_df = read.csv(probeMap_path,sep='\t')
###
probeMap_df <- probeMap_df %>%  dplyr::select(1,2)
colnames(probeMap_df) <- c("Ensembl_ID","gene_name")

exp_mat <- inner_join(probeMap_df, exp_mat, by = "Ensembl_ID") %>% dplyr::select(-1)
exp_mat <- aggregate(.~ gene_name, exp_mat, mean)
exp_mat <- exp_mat %>% column_to_rownames("gene_name")

# =========== TCGA deseq2 normalization


deseq2_exp_mat = exp_mat
condition_table = factor(c(rep('sen',608),rep('unsen',609)))
# deseq2_exp_mat = as.matrix(exp_mat%>%column_to_rownames('gene_id'))
deseq2_exp_mat = apply(deseq2_exp_mat,2,as.integer)
rownames(deseq2_exp_mat) = rownames(exp_mat)

deseq2_normalized_exp_mat = deseq2_normalization(deseq2_exp_mat,condition_table)


# === select pathway gene
exp_mat_cc = deseq2_normalized_exp_mat[unique(pw_list_df$F),] # 选择pathway基因
exp_mat_cc = na.omit(exp_mat_cc)


###############
# 3 pathway PCA
###############

prot_list = list()
for(item in unique(pw_list_df$A)){
    prot_list[[item]] = unique(as.vector((pw_list_df%>%filter(A==item))$F))
}
pw_g_list = prot_list
exp_mat_gsva = deseq2_normalized_exp_mat[unique(pw_list_df$F),] # 选择pathway基因
exp_mat_gsva = na.omit(exp_mat_gsva)
gsva_matrix <- gsva(as.matrix(exp_mat_gsva), pw_g_list, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE, parallel.sz=1)



library(rgl)
library(plotly)
attach(mtcars)

df_for_plot = as.data.frame(t(gsva_matrix))

plot3d(df_for_plot$UPS,df_for_plot$CMA,df_for_plot$UPR,clo="red",size=5)
dir.create('result2')
htmlwidgets::saveWidget(rglwidget(width = 520, height = 520), 
                        file = "./result2/3dscatter.html",
                        libdir = "libs",
                        selfcontained = FALSE
                        )

plot_ly(data=df_for_plot, x = ~UPS, y = ~CMA, z = ~UPR, size=1)


# 研究一个函数进行全部比较的方法
cor.test(df_for_plot$UPS, df_for_plot$UPR)
cor.test(df_for_plot$UPS, df_for_plot$CMA)
cor.test(df_for_plot$CMA, df_for_plot$UPR)


#####################
# census clustering
#####################

# =================

exp_mat_cc = deseq2_normalized_exp_mat[unique(pw_list_df$F),] # 选择pathway基因
exp_mat_cc = na.omit(exp_mat_cc)

# 对样本聚类
# Consensus Cluster #################
# input of ConsensusClusterPlus: gene×sample ,this func is do the cluster for sample
##
results = ConsensusClusterPlus(t(scale(t(exp_mat_cc))),maxK=8)
save(results, file = sprintf( paste(c(result_dir,sprintf('consenCluster_res_%s_for_sample.RData',project_name)), collapse='/') ))

## post cluster analysis

# pw_sort_ls = read.table('./pathway_sort_ls.txt',header = FALSE,sep = '\n')

### 这里是取4类的结果来分析
class_df = as.data.frame(results[[7]][['consensusClass']],col.names=c('class'))
colnames(class_df)=c('class')
class_df['sample'] = rownames(class_df)

tmp1 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['1']]
tmp2 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['2']]
tmp3 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['3']]
tmp4 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['4']]
tmp5 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['5']]
tmp6 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['6']]
tmp7 = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['7']]

annotation_col = data.frame(
  Cluster_col = factor(c(rep('c1',length(tmp1)) , rep('c2',length(tmp2)) , rep('c3',length(tmp3)) , rep('c4',length(tmp4)), rep('c5',length(tmp5)), rep('c6',length(tmp6)),
                    rep('c7',length(tmp7))))
)
rownames(annotation_col) = c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7)


# 对基因聚类
# Consensus Cluster #################
# input of ConsensusClusterPlus: gene×sample ,this func is do the cluster for sample
##
results_for_gene = ConsensusClusterPlus(scale(t(exp_mat_cc)),maxK=8)
save(results_for_gene, file = sprintf( paste(c(result_dir,sprintf('consenCluster_res_%s_for_gene.RData',project_name)), collapse='/') ))

## post cluster analysis

# pw_sort_ls = read.table('./pathway_sort_ls.txt',header = FALSE,sep = '\n')

### 这里是取6类的结果来分析
class_df = as.data.frame(results_for_gene[[6]][['consensusClass']],col.names=c('class'))
colnames(class_df)=c('class')
class_df['sample'] = rownames(class_df)

tmp1_g = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['1']]
tmp2_g = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['2']]
tmp3_g = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['3']]
tmp4_g = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['4']]
tmp5_g = split(as.matrix(class_df)[,2], as.matrix(class_df[,1]))[['5']]

annotation_row = data.frame(
  Cluster_row = factor(c(rep('c1_g',length(tmp1_g)) , rep('c2_g',length(tmp2_g)) , rep('c3_g',length(tmp3_g)) , rep('c4_g',length(tmp4_g)), rep('c5_g',length(tmp5_g))))
)
rownames(annotation_row) = c(tmp1_g,tmp2_g,tmp3_g,tmp4_g,tmp5_g)


# ======= enrichment
library("clusterProfiler")
gmtfile <- '~/my_project/data/pathway/c2.cp.kegg.v7.2.symbols.gmt'
c5 <- read.gmt(gmtfile);c5

# low level enrichment
tmp_enrich_df = pw_list_df[,c('B','F')]
colnames(tmp_enrich_df) = c('term','gene')

kk <- enricher(as.vector(tmp1_g), TERM2GENE=tmp_enrich_df, 
                    pAdjustMethod = "BH", pvalueCutoff = 1, 
                    minGSSize = 0, maxGSSize = 500)
as.data.frame(kk)

# top level enrichment
tmp_enrich_df = pw_list_df[,c('A','F')]
colnames(tmp_enrich_df) = c('term','gene')

kk <- enricher(as.vector(tmp1_g), TERM2GENE=tmp_enrich_df, 
                    pAdjustMethod = "BH", pvalueCutoff = 1, 
                    minGSSize = 0, maxGSSize = 500)
as.data.frame(kk)



# ======= Visualization

# heatmap

## 对通路的注释的整理
anno_ls = c()
for(g in c(tmp1_g,tmp2_g,tmp3_g,tmp4_g,tmp5_g)){
    tmp = paste(unique((pw_list_df%>%filter(F==g))$A),collapse='_')
    if(tmp==''){
        anno_ls = c(anno_ls, 'X')
    }else{
        anno_ls = c(anno_ls, tmp)
    }
}

annotation_row = data.frame(
  Cluster_row = factor(c(rep('c1_g',length(tmp1_g)) , rep('c2_g',length(tmp2_g)) , rep('c3_g',length(tmp3_g)) , rep('c4_g',length(tmp4_g)), rep('c5_g',length(tmp5_g)))),
  term_row = factor(anno_ls)
)
rownames(annotation_row) = c(tmp1_g,tmp2_g,tmp3_g,tmp4_g,tmp5_g)

annotation_col = data.frame(
  Cluster_col = factor(c(rep('c1',length(tmp1)) , rep('c2',length(tmp2)) , rep('c3',length(tmp3)) , rep('c4',length(tmp4)), rep('c5',length(tmp5)), rep('c6',length(tmp6)),
                    rep('c7',length(tmp7))))
)
rownames(annotation_col) = c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7)

# 控制柱饰条颜色
anno_colors = list(
            Cluster_col = c(c1="red",c2="blue",c3="green",c4="pink",c5='black',c6='yellow',c7='grey'),
            Cluster_row = c(c1_g="red",c2_g="blue",c3_g="green",c4_g="pink",c5_g='black'),
            term_row = c(CMA='Cyan',UPR='Gold',UPS='Coral',UPS_CMA='DimGray',X='DarkViolet')
            )


# 控制极端值，保证热图好看
exp_mat_cc_tmp = t(scale(t(exp_mat_cc)))
exp_mat_cc_tmp[exp_mat_cc_tmp>5] = 5
exp_mat_cc_tmp[exp_mat_cc_tmp<-5] = -5

# 图例色条调整
bk <- c(seq(-5,-0.1,by=0.01),seq(0,5,by=0.01))


pheatmap(
    exp_mat_cc_tmp[c(tmp1_g,tmp2_g,tmp3_g,tmp4_g,tmp5_g),
               c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7)],
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    #scale = 'row',
    #color = colorRampPalette(c("blue","white", "red"))(50),
    color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
    legend_breaks=seq(-5,5,2.5),
    breaks=bk,
    # filename='./cencluster_hmap_ccgene_ccsample_modi_color.pdf',
    filename = paste(c(result_dir,sprintf('cencluster_hmap_ccgene_ccsample_%s.pdf',project_name)), collapse='/')
    cellwidth = 2, cellheight = 10,
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = anno_colors,
    border_color = "grey"
    )























