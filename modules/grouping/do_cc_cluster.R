


pn_top = 0.15

gene_var = apply(normalized_exp_mat%>%column_to_rownames('gene_symbol'),1,var)
sel_genes = names(gene_var[base::order(gene_var,decreasing=TRUE)][1:round(length(gene_var)*pn_top)])

exp_mat_cc = normalized_exp_mat%>%filter(gene_symbol%in%sel_genes)%>%column_to_rownames('gene_symbol')


# 对样本聚类
# Consensus Cluster #################
# input of ConsensusClusterPlus: gene×sample ,this func is do the cluster for sample
##
results = ConsensusClusterPlus(t(scale(t(exp_mat_cc))),maxK=8,clusterAlg='km',pItem=0.8,reps=50)
#save(results, file = sprintf( paste(c(result_dir,sprintf('consenCluster_res_%s_for_sample.RData',paste(project_name,pn_top,sep='_'))), collapse='/') ))




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














