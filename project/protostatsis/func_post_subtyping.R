# group_df$sample_id_B 是针对第二组学的样本名与第一组学（rnaseq）样本名不一样的情况的，提供第二组学样本名来对应
# compare_drug_sen_two_drug检验靶点的富集情况！！！
#


PATH_SRC='/home/lgh/my_project/omic/'
PATH_DATA='/home/lgh/my_project/data/'
source(sprintf('%s/modules/clean_data/normalization.R',PATH_SRC))
source(sprintf('%s/modules/bioinformatics_wrapper/enrich.R',PATH_SRC))
source(sprintf('%s/modules/bioinformatics_wrapper/EDA.R',PATH_SRC))
source('/home/lgh/my_project/omic/modules/bioinformatics_wrapper/methy.R')

#  group_list = list('g2'=tmp2,'g3'=tmp3,'g5'=tmp5)
#  comb_res = combn(names(group_list),2)
#> g2	g2	g3
#> g3	g5	g5

# group_df: sample;group
#
#

##############################################
# do eda by group and GSEA
##############################################
# comp
## compare
## samples
## deseq2_res
## gsea_df_go/gsea_df_kegg/gsea_df_react

# exp_mat: deseq2给count，limma给norm_exp_mat
# ~~~group_list  group_df

compare_transcriptome <- function(exp_mat,gmtfile,group_df,method='deseq2',gsea=FALSE){
  cat(sprintf('use %s method!\n\n',method))

  group_comp_res = list()
  comb_res = combn(unique(as.vector(group_df$group)),2)

  # i = 1 # for(i in 1:dim(comb_res)[2])
  for(i in 1:dim(comb_res)[2]){
      
    g1 = comb_res[,i][1]
    g2 = comb_res[,i][2]
    print(paste(c('# # now is',i,'compare between',g1,'&',g2,'...'),collapse=' '))

    comp = paste('comp',g1,g2,sep='_')
    group_comp_res[[comp]] = list()

    group_comp_res[[comp]]$compare = c(g1,g2)

    g1_samples = (group_df%>%filter(group==g1))$sample
    g2_samples = (group_df%>%filter(group==g2))$sample
    group_comp_res[[comp]]$samples = list(g1_samples=g1_samples,g2_samples=g2_samples)

    g1_num = length(g1_samples)
    g2_num = length(g2_samples)
    # return(exp_mat)
    # deseq2 =====================================
    exp_mat_for_dea = exp_mat[,c(g1_samples,g2_samples)]
    condition_table = c(rep(g1, g1_num), rep(g2, g2_num))
    # ===

    ## do DEA
    ## Row.names,stat,padj,log2FoldChange
    if(method=='deseq2'){
      rn = rownames(exp_mat_for_dea)
      exp_mat_for_dea = apply(exp_mat_for_dea,2,as.integer)
      rownames(exp_mat_for_dea) = rn
      group_comp_res[[comp]]$dea_res = do_deseq2(exp_mat_for_dea, condition_table)$resdata
    }else{
      group_comp_res[[comp]]$dea_res = do_limma(exp_mat_for_dea, condition_table)$resdata
    }


    signif_dea_res = group_comp_res[[comp]]$dea_res%>%filter(padj<0.01 & log2FoldChange>2) ; print(paste('signif deseq2 res ',dim(signif_dea_res)[1], sep='') ) ###

    # do gsea  =====================================
    group_comp_res[[comp]]$gsea_res = list()
    if(gsea==TRUE){
      dea_res = group_comp_res[[comp]]$dea_res
      #return(dea_res)
      for(gsea_type in names(gmtfile)){
        gene_df = data.frame(gene=dea_res$row.names, score_for_rank=dea_res$stat)  ###
        #gmtfile = paste(PATH_DATA,'/pathway/c5.go.bp.v7.2.symbols.gmt',sep=''  )            ### go all
        group_comp_res[[comp]][['gsea_res']][[gsea_type]] = GSEA_enrich(gene_df, gmtfile[[gsea_type]], pvalueCutoff=1)      
      }

      #gene_df = data.frame(gene=dea_res$row.names, score_for_rank=dea_res$stat)  ###
      #gmtfile = paste(PATH_DATA,'/pathway/c5.go.bp.v7.2.symbols.gmt',sep=''  )            ### go all
      #group_comp_res[[comp]]$gsea_df_go = GSEA_enrich(gene_df, gmtfile$go, pvalueCutoff=0.1)

      #gene_df = data.frame(gene=dea_res$row.names, score_for_rank=dea_res$stat)  ###
      #gmtfile = paste(PATH_DATA,'/pathway/c2.cp.kegg.v7.2.symbols.gmt',sep='')            ### kegg
      #group_comp_res[[comp]]$gsea_df_kegg = GSEA_enrich(gene_df, gmtfile$kegg, pvalueCutoff=0.1)

      #gene_df = data.frame(gene=dea_res$row.names, score_for_rank=dea_res$stat)  ###
      #gmtfile = paste(PATH_DATA,'/pathway/reactome.gmt',sep='')                           ### reactome
      #group_comp_res[[comp]]$gsea_df_react = GSEA_enrich(gene_df, gmtfile$reactome, pvalueCutoff=0.1)
    }

  }
      
  return(group_comp_res)
}

if(FALSE){
  gmtfile = list(
      go=paste(PATH_DATA,'/pathway/c5.go.bp.v7.2.symbols.gmt',sep=''  ),
      kegg=paste(PATH_DATA,'/pathway/c2.cp.kegg.v7.2.symbols.gmt',sep=''),
      reactome=paste(PATH_DATA,'/pathway/reactome.gmt',sep='')
    )
  compare_transcriptome_res = compare_transcriptome(exp_mat=normalized_exp_mat,
                        gmtfile=gmtfile,
                        group_df=subtyping_res$sample_group_df,method='limma',gsea=FALSE)
}
## 



##############################################
# survival
##############################################

## 有没有多组生存分析的方法

# surv_df = data.frame(OS=,OSS=,group=)  ###
## 前面创建一个有生存和分组的数据框

if(FALSE){ # TCGA
  group_df = data.frame(sample=c(group_list[['g2']],group_list[['g3']],group_list[['g5']]),
                        group=c(rep('g2',length(group_list[['g2']])),rep('g3',length(group_list[['g3']])),
                                rep('g5',length(group_list[['g5']]))
                              ) 
                        )
  surv_df = read_tsv('~/my_project/data/TCGA/BRCA/TCGA-BRCA.survival.tsv')
  colnames(surv_df) = c('sample','OSS','_PATIENT','OS')
  # source('/home/lgh/my_project/omic/modules/plot/survival.R')
}
if(FALSE){ # VIZOME
  cli_df = read_tsv('~/my_project/data/VIZOME/nature_aml_clinical_annotation.txt')
  cli_df[['LabId...1']] = paste('X',str_replace(cli_df[['LabId...1']],'-','.'),sep='')

  tmp_OSS = cli_df$causeOfDeath
  tmp_OSS[tmp_OSS%in%c('Dead-Disease','Dead-Treatment')] = 1
  tmp_OSS[tmp_OSS%in%c('Dead-Unknown','Dead-Other','Alive')] = 0
  suv_df = data.frame(sample=cli_df$LabId...1,
                      OS=cli_df$overallSurvival, 
                      OSS=as.numeric(tmp_OSS) )

  compare_survival_res = compare_survival(suv_df, group_df=sample_group_df,plot=FALSE)
  compare_survival_res$plot

  compare_survival_res = compare_survival(suv_df, group_df=sample_group_df%>%filter(group%in%c('c1','c2','c3')),plot=FALSE)
  compare_survival_res$plot
}
##

# 要获得生存分析的关键
#
compare_survival <- function(surv_df, group_df,plot=TRUE){
  source('/home/lgh/my_project/omic/modules/plot/survival.R')

  df_for_suv_plot = merge(surv_df,group_df,by='sample')
  var = 'group'
  res=plot_survival_binary(df_for_suv_plot, var, plot=TRUE)
  if(plot){
    plot(res$plot)
  }
  return(res)
}

#tmp_surv_df = surv_df_all%>%filter(sample%in%c(g1_samples, g2_samples))
#var = 'group'
#surv_diff <- survdiff(Surv(tmp_surv_df[['OS']], tmp_surv_df[['OSS']]) ~ tmp_surv_df[[var]])
#pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)   




##############################################
# 免疫浸润分析
##############################################
if(FALSE){
  infil_df = read_csv('~/my_project/data/TCGA/infiltration_estimation_for_tcga.csv.gz')
  group_df = data.frame(sample=c(group_list[['g2']],group_list[['g3']],group_list[['g5']]),
                        group=c(rep('g2',length(group_list[['g2']])),rep('g3',length(group_list[['g3']])),
                                rep('g5',length(group_list[['g5']]))
                              ) 
                        )
  group_df$sample_id_B = sapply(group_df$sample,function(x) substr(x,1,nchar(x)-1))
  # comb_res
}
##
##
compare_immu <- function(infil_df,group_df,comb_res){
  df_for_comp = merge(group_df,
                  reshape2::melt(infil_df, id.vars = c("cell_type")),
                  by.x='sample_id_B',by.y='cell_type')

  df_for_comp$immu_cell = sapply(strsplit(as.vector(df_for_comp$variable),'_'),'[[',1)
  df_for_comp$immu_algo = sapply(strsplit(as.vector(df_for_comp$variable),'_'),'[[',2)

  # cal the signif immu
  comp_res_df = compare_means(formula=value~group,data=df_for_comp,group.by='variable',p.adjust.method='fdr')

  comp_res_df$immu_cell = sapply(strsplit(as.vector(comp_res_df$variable),'_'),'[[',1)
  comp_res_df$immu_algo = sapply(strsplit(as.vector(comp_res_df$variable),'_'),'[[',2)

  # filter
  signif_comp_res_df = comp_res_df%>%filter(p.adj<0.01)

  # get infil obj present in 5 algo
  signif_immu_cell = c()
  for(c in unique(signif_comp_res_df$immu_cell)){
      #print(length((signif_comp_res_df%>%filter(immu_cell==c))$immu_algo))
      if( length(unique((signif_comp_res_df%>%filter(immu_cell==c))$immu_algo))>=5 ){
          signif_immu_cell = c(signif_immu_cell,c)
      }
  }

  # plot
  my_comparisons = split(comb_res, rep(1:ncol(comb_res), each = nrow(comb_res)))

  p <- ggboxplot(df_for_comp%>%filter(immu_cell%in%signif_immu_cell), 
            x = "group", y = "value",
            color = "group", palette = "nature", facet.by = "variable")+
            theme(axis.text.x = element_text(angle = 75, vjust=0.5)) + stat_compare_means(aes(group=group),comparisons=my_comparisons)
            #stat_compare_means(aes(group=variable),comparisons=my_comparisons)
  plot(p)

  return(list(df_for_comp=df_for_comp,
              comp_res_df=comp_res_df,
              signif_immu_cell=signif_immu_cell,
              plot=p))
}

##############################################
# 多组学关联分析
## 关联甲基化
##############################################

compare_methy <- function(methy_dat){

  library(data.table)
  library(impute)
  library(ChAMP)
  library(stringr)
  library(tibble)

  myNorm = methy_dat$myNorm
  pd = methy_dat$pd

  ## 探索性分析 pca
  library(FactoMineR)
  library(factoextra) 

  dat <- t(myNorm)

  table(pd$group)

  dat.pca <- PCA(dat, graph = FALSE) 
  fviz_pca_ind(dat.pca,
              geom.ind = "point", 
              col.ind = pd$group, 
              addEllipses = TRUE, 
              legend.title = "Groups"
              )

  ## 探索性分析 聚类热图
  library(pheatmap)
  cg=names(tail(sort(apply(myNorm,1,sd)),1000))
  library(pheatmap)
  ac=data.frame(group=pd$group)
  rownames(ac)=colnames(myNorm)  
  pheatmap(myNorm[cg,],show_colnames =F,show_rownames = F,
          annotation_col=ac)

  ## 探索性分析 相关性矩阵聚类
  pheatmap::pheatmap(cor(myNorm[cg,]),
                    annotation_col = ac,
                    annotation_row = ac,
                    show_rownames = F,
                    show_colnames = F)

  ## 组间差异甲基化

  #group_list <- pd$group
  #group_list[is.na(group_list)] = 'NN'

  myDMP <- champ.DMP(beta = myNorm[,pd$group!='N'],pheno=as.vector(pd$group[pd$group!='N']))

  return(myDMP)
}

if(FALSE){
  group_df = data.frame(sample=c(group_list[['g2']],group_list[['g3']],group_list[['g5']]),
                        group=c(rep('g2',length(group_list[['g2']])),rep('g3',length(group_list[['g3']])),
                                rep('g5',length(group_list[['g5']]))
                              ) 
                        )
  group_df$sample_id_B = sapply(group_df$sample,function(x) substr(x,1,nchar(x)-1))

  methy_df = data.table::fread("~/my_project/data/TCGA/BRCA/TCGA-BRCA.methylation450.tsv.gz", data.table = F)
  methy_df = column_to_rownames(methy_df,"Composite Element REF")
  colnames(methy_df)= str_sub(colnames(methy_df),1,15)

  methy_dat = load_TCGA_methy_and_preprocess(group_df, methy_path)
}
##
##

### 差异可视化
if(FALSE){
  tmp_df = cbind(data.frame(myNorm['cg12098441',]),pd)
  colnames(tmp_df)[1] = 'value'
  tmp_df = tmp_df%>%filter(group!='-1')
  ggboxplot(tmp_df,x='group',y='value',add='jitter',color='group') + stat_compare_means()
}

# return(myDMP)

##############################################
# 多组学关联分析
## miRNA  ### 多组互相比较问题
##############################################
if(FALSE){
  mirna_mat = read_tsv('~/my_project/data/TCGA/BRCA/TCGA-BRCA.mirna.tsv.gz')
  group_df = data.frame(sample=c(group_list[['g2']],group_list[['g3']],group_list[['g5']]),
                        group=c(rep('g2',length(group_list[['g2']])),rep('g3',length(group_list[['g3']])),
                                rep('g5',length(group_list[['g5']]))
                              ) 
                        )
  group_df$sample_id_B = group_df$sample # miRNA的sample name与rnaseq一致
}
##
##
compare_miRNA <- function(mirna_mat, group_df){
  library(limma)

  mirna_mat = mirna_mat[,]%>%column_to_rownames('miRNA_ID')
  tmp = group_df%>%filter(sample_id_B%in%colnames(mirna_mat))
  mat_for_limma = mirna_mat[,tmp$sample]
  condition_table = tmp$group

  limma_res = do_limma(mat_for_limma, condition_table)

  ggscatter(limma_res,x='logFC',y='adj.P.Val') # 火山图展示

  return(limma_res)
}

##############################################
# 多组学关联分析
## 突变
##############################################

load_TCGA_mut_2_maf <- function(maf.dir){
  if('data.frame'%in%class(maf.dir)){
    maf.df = maf.dir
  }else{
    maf.df = read_tsv(maf.dir)
  }
  ##
  ## preprocess snv file, make maf
  colnames(maf.df) =c( "Tumor_Sample_Barcode", "Hugo_Symbol", 
                    "Chromosome", "Start_Position", 
                  "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
                  "amino_acid_change" , 'effect' ,"Consequence",
                  "vaf" ) 
  maf.df$Entrez_Gene_Id =1
  maf.df$Center ='ucsc'
  maf.df$NCBI_Build ='GRCh38'
  maf.df$NCBI_Build ='GRCh38'
  maf.df$Strand ='+'
  maf.df$Variant_Classification = maf.df$effect
  tail(sort(table(maf.df$Variant_Classification )))
  maf.df$Tumor_Seq_Allele1 = maf.df$Reference_Allele
  maf.df$Variant_Type = ifelse(
    maf.df$Reference_Allele %in% c('A','C','T','G') & maf.df$Tumor_Seq_Allele2 %in% c('A','C','T','G'),
    'SNP','INDEL'
  )
  table(maf.df$Variant_Type )
  #maf = read.maf(maf = maf.df,
  #                    vc_nonSyn=names(tail(sort(table(maf.df$Variant_Classification )))))
  return(maf.df)
}

pair_compare_maf <- function(maf.df, group_df, m1Name, m2Name, m1show, m2show, my_genes=NULL, sig_th=0.05, plot=FALSE){
  library(maftools)
  colnames(group_df)[1] = 'Tumor_Sample_Barcode'
  maf.df1=maf.df%>%filter(Tumor_Sample_Barcode%in%(group_df%>%filter(group==m1Name))$Tumor_Sample_Barcode)
  maf1 = read.maf(maf = maf.df1,
                      vc_nonSyn=names(tail(sort(table(maf.df1$Variant_Classification )))),verbose = FALSE)
  maf.df2=maf.df%>%filter(Tumor_Sample_Barcode%in%(group_df%>%filter(group==m2Name))$Tumor_Sample_Barcode)
  maf2 = read.maf(maf = maf.df2,
                      vc_nonSyn=names(tail(sort(table(maf.df2$Variant_Classification )))),verbose = FALSE)

  ### 寻找显著差异突变的基因
  pt.vs.rt <- mafCompare(m1 = maf1, m2 = maf2, m1Name = m1Name, m2Name = m2Name, minMut = 1)
  tmb_df = rbind(data.frame(tmb(maf = maf1), group=m1Name),data.frame(tmb(maf = maf2), group=m2Name))
  if(plot){
    # print(pt.vs.rt)
    ### 森林图展示显著突变基因
    forestPlot(mafCompareRes = pt.vs.rt, pVal = sig_th, color = c('royalblue', 'maroon'), geneFontSize = 0.8)
    ### 瀑布图比较两组的突变情况
    genes = pt.vs.rt$result$Hugo_Symbol[1:20]
    if(!is.null(my_genes)){
      genes = c(my_genes,genes)
    }
    coOncoplot(m1 = maf1, m2 = maf2, m1Name = m1show, m2Name = m2show, genes = genes, removeNonMutated = TRUE)
    ## lollipopPlot2图展示某基因上面两组样本的突变情况
    #lollipopPlot2(m1 = brca.maf.g2, m2 = brca.maf.g3, gene = "PIK3CA", 
    #              AACol1 = "amino_acid_change", AACol2 = "amino_acid_change", m1_name = "g1", m2_name = "g2")
    ### 展示两组的突变通路情况
    #OncogenicPathways(maf = brca.maf.g2)
    #OncogenicPathways(maf = brca.maf.g3)
    ### 两组计算driver gene
    #brca.maf.g2.sig = oncodrive(maf = brca.maf.g2, AACol = 'amino_acid_change', minMut = 5, pvalMethod = 'zscore')
    #head(brca.maf.g2.sig)
    #plotOncodrive(res = brca.maf.g2.sig, fdrCutOff = 0.05, useFraction = TRUE)
    ### 两组样本的drugInteractions
    #dgi = drugInteractions(maf = brca.maf.g2, fontSize = 0.75)
    #dgi = drugInteractions(maf = brca.maf.g3, fontSize = 0.75)
    ### 两组的突变相关性
    somaticInteractions(maf = maf1, top = 25, pvalue = c(0.01, 0.05))
    somaticInteractions(maf = maf2, top = 25, pvalue = c(0.01, 0.05))
    ### 比较两组的tmb
    library(ggpubr)
    # tmb_df = rbind(data.frame(tmb(maf = maf1), group=m1Name),data.frame(tmb(maf = maf2), group=m2Name))
    p=ggboxplot(tmb_df,x='group',y='total_perMB_log',color='group')+stat_compare_means()
    plot(p)
  }
  return(list(pt.vs.rt=pt.vs.rt, tmb_df=tmb_df, maf1=maf1, maf2=maf2))

}

compare_maf <- function(maf.df,group_df){
  library(maftools)
  colnames(group_df)[1] = 'Tumor_Sample_Barcode'
  maf = read.maf(maf = maf.df, clinicalData = group_df[,c(1,2)],
                vc_nonSyn=names(tail(sort(table(maf.df$Variant_Classification )))),verbose = FALSE)
  oncoplot(
    maf = maf,
    clinicalFeatures = 'group',
    sortByAnnotation = TRUE)

  fab.ce = clinicalEnrichment(maf = maf, clinicalFeature = 'group')

  return(list(fab.ce=fab.ce,
            maf=maf))
}

pair_compare_maf_for_depmap_cellline_name <- function(info_df, maf.df, group_df, m1Name, m2Name, 
                                                m1show, m2show, sig_th=0.05, plot=FALSE){ # 待完善

  group_df = merge(info_df[,c('DepMap_ID','stripped_cell_line_name')],
                group_df,by.x='DepMap_ID',by.y='Tumor_Sample_Barcode')[-1]
  colnames(group_df)[1] = 'Tumor_Sample_Barcode'

  maf.df = merge(info_df[,c('DepMap_ID','stripped_cell_line_name')],maf.df,by.x='DepMap_ID',by.y='DepMap_ID')[-1]
  maf.df$Tumor_Sample_Barcode = maf.df$stripped_cell_line_name


  library(maftools)
  colnames(group_df)[1] = 'Tumor_Sample_Barcode'
  maf.df1=maf.df%>%filter(Tumor_Sample_Barcode%in%(group_df%>%filter(group==m1Name))$Tumor_Sample_Barcode)
  maf1 = read.maf(maf = maf.df1,
                      vc_nonSyn=names(tail(sort(table(maf.df1$Variant_Classification )))),verbose = FALSE)
  maf.df2=maf.df%>%filter(Tumor_Sample_Barcode%in%(group_df%>%filter(group==m2Name))$Tumor_Sample_Barcode)
  maf2 = read.maf(maf = maf.df2,
                      vc_nonSyn=names(tail(sort(table(maf.df2$Variant_Classification )))),verbose = FALSE)

  ### 寻找显著差异突变的基因
  pt.vs.rt <- mafCompare(m1 = maf1, m2 = maf2, m1Name = m1Name, m2Name = m2Name, minMut = 1)
  tmb_df = rbind(data.frame(tmb(maf = maf1), group=m1Name),data.frame(tmb(maf = maf2), group=m2Name))
  if(plot){
    # print(pt.vs.rt)
    ### 森林图展示显著突变基因
    forestPlot(mafCompareRes = pt.vs.rt, pVal = sig_th, color = c('royalblue', 'maroon'), geneFontSize = 0.8)
    ### 瀑布图比较两组的突变情况
    genes = pt.vs.rt$result$Hugo_Symbol[1:20]
    coOncoplot(m1 = maf1, m2 = maf2, m1Name = m1show, m2Name = m2show, genes = genes, removeNonMutated = TRUE)
    ## lollipopPlot2图展示某基因上面两组样本的突变情况
    #lollipopPlot2(m1 = brca.maf.g2, m2 = brca.maf.g3, gene = "PIK3CA", 
    #              AACol1 = "amino_acid_change", AACol2 = "amino_acid_change", m1_name = "g1", m2_name = "g2")
    ### 展示两组的突变通路情况
    #OncogenicPathways(maf = brca.maf.g2)
    #OncogenicPathways(maf = brca.maf.g3)
    ### 两组计算driver gene
    #brca.maf.g2.sig = oncodrive(maf = brca.maf.g2, AACol = 'amino_acid_change', minMut = 5, pvalMethod = 'zscore')
    #head(brca.maf.g2.sig)
    #plotOncodrive(res = brca.maf.g2.sig, fdrCutOff = 0.05, useFraction = TRUE)
    ### 两组样本的drugInteractions
    #dgi = drugInteractions(maf = brca.maf.g2, fontSize = 0.75)
    #dgi = drugInteractions(maf = brca.maf.g3, fontSize = 0.75)
    ### 两组的突变相关性
    somaticInteractions(maf = maf1, top = 25, pvalue = c(0.01, 0.05))
    somaticInteractions(maf = maf2, top = 25, pvalue = c(0.01, 0.05))
    ### 比较两组的tmb
    library(ggpubr)
    # tmb_df = rbind(data.frame(tmb(maf = maf1), group=m1Name),data.frame(tmb(maf = maf2), group=m2Name))
    p=ggboxplot(tmb_df,x='group',y='total_perMB_log',color='group')+stat_compare_means()
    plot(p)
  }
  return(list(pt.vs.rt=pt.vs.rt, tmb_df=tmb_df, maf1=maf1, maf2=maf2))
  # 出来后可以个性化绘图
  # pdf('./mafplot.pdf',width=10,height=10)
  # genes=pair_compare_res$pt.vs.rt$result$Hugo_Symbol[1:10]
  # coOncoplot(m1 = pair_compare_res$maf1, m2 = pair_compare_res$maf2, m1Name = 'sen', m2Name = 'unsen', genes = genes, removeNonMutated = TRUE,showSampleNames=TRUE,SampleNamefont=1)
  # dev.off()
}

if(FALSE){
  library(maftools)
  group_df = data.frame(sample=c(group_list[['g2']],group_list[['g3']],group_list[['g5']]),
                        group=c(rep('g2',length(group_list[['g2']])),rep('g3',length(group_list[['g3']])),
                                rep('g5',length(group_list[['g5']]))
                              ) 
                        )
  group_df$sample_id_B = sapply(group_df$sample,function(x) substr(x,1,nchar(x)-1))

  maf.dir = '/home/lgh/my_project/data/TCGA/BRCA/TCGA-BRCA.mutect2_snv.tsv.gz'
}


if(FALSE){
  brca.maf.df = load_TCGA_mut_2_maf(maf.dir)

  ###
  ###
  maf = read.maf(maf = maf.df,
                vc_nonSyn=names(tail(sort(table(maf.df$Variant_Classification )))))
  oncoplot(maf = brca.maf.df,top = 10) # 高频突变的前10个基因

  ## g2 g3之间进行比较

  compare_maf(maf.df=brca.maf.df, group_df=group_df, m1Name='g2', m2Name='g3')
}

if(FALSE){ # 合起来oncoplot和多重比较
  
  colnames(group_df)[1] = 'Tumor_Sample_Barcode'
  maf = read.maf(maf = maf.df, clinicalData = group_df[,c(1,2)],
                vc_nonSyn=names(tail(sort(table(maf.df$Variant_Classification )))))
  oncoplot(
    maf = maf,
    clinicalFeatures = 'group',
    sortByAnnotation = TRUE)

  fab.ce = clinicalEnrichment(maf = maf, clinicalFeature = 'group')
  plotEnrichmentResults(enrich_res = fab.ce, 
    pVal = 0.05, 
    geneFontSize = 0.7, 
    annoFontSize = 0.7,
    legendFontSize = 1)

  #oncostrip(maf, genes = NULL, sort = TRUE, annotation = NULL,
  #  removeNonMutated = FALSE, top = 5, showTumorSampleBarcodes = FALSE,
  #  colors = NULL)

}


##############################################
# 药敏
##############################################
if(FALSE){ # old version
  drug_family_df = read_tsv('~/my_project/data/VIZOME/nature_aml_drug_families.txt')
  drug_anno_df = read_tsv('~/my_project/data/VIZOME/nature_aml_drug_annotation.txt')
  drug_sen_df_viz = read_tsv('~/my_project/data/VIZOME/nature_aml_drug_sen.tsv')
  ###
  df_for_comp_viz = merge(sample_group_df,drug_sen_df_viz,by.x='sample',by.y='lab_id')
  comp_res_ic_viz = compare_means(ic50~group,df_for_comp_viz,group.by='inhibitor',p.adjust.method='fdr')
  comp_res_ic_viz[comp_res_ic_viz$p.adj<0.05,]

  cal_mean = df_for_comp_viz%>%
    group_by(group,inhibitor)%>%
    summarise(y=mean(ic50))
  colnames(cal_mean) = c('group1','inhibitor','mean_g1')
  comp_res_ic_viz = merge(comp_res_ic_viz,cal_mean, by=c('inhibitor','group1'))
  colnames(cal_mean) = c('group2','inhibitor','mean_g2')
  comp_res_ic_viz = merge(comp_res_ic_viz,cal_mean, by=c('inhibitor','group2'))

  comp_res_ic_viz = merge(comp_res_ic_viz,drug_anno_df,by='inhibitor') ## 不是每个药都有注释

}

# drug_sen_df
# group_df
# drug_anno_df
# ic50auc
# 
compare_drug_sen <- function(drug_sen_df, group_df, drug_anno_df=NULL, ic50auc='i50'){
  df_for_comp = merge(group_df,drug_sen_df,by.x='sample',by.y='sample')
  comp_res = compare_means(ic50~group,df_for_comp,group.by='inhibitor',p.adjust.method='fdr')
  comp_res = comp_res[order(comp_res$p.adj),]

  if(ic50auc=='ic50'){

    comp_res = compare_means(ic50~group,df_for_comp,group.by='inhibitor',p.adjust.method='fdr')
    comp_res = comp_res[order(comp_res$p.adj),]

    cal_mean = df_for_comp%>%
      group_by(group,inhibitor)%>%
      summarise(y=mean(ic50,na.rm = TRUE))
  }else{

    comp_res = compare_means(auc~group,df_for_comp,group.by='inhibitor',p.adjust.method='fdr')
    comp_res = comp_res[order(comp_res$p.adj),]

    cal_mean = df_for_comp%>%
      group_by(group,inhibitor)%>%
      summarise(y=mean(auc,na.rm = TRUE))
  }

  colnames(cal_mean) = c('group2','inhibitor','mean_g2')
  comp_res = merge(comp_res,cal_mean, by=c('inhibitor','group2'))
  colnames(cal_mean) = c('group1','inhibitor','mean_g1')
  comp_res = merge(comp_res,cal_mean, by=c('inhibitor','group1'))
  comp_res=comp_res[,c(c(1:9),11,10)]

  comp_res = comp_res[order(comp_res$p),]

  comp_res_anno = NULL
  if(!is.null(drug_anno_df)){
    comp_res_anno = merge(comp_res,
                  drug_anno_df,by='inhibitor',
                  all.x=TRUE) ## 不是每个药都有注释,但每个药都要保留
  }
  return(list(comp_res=comp_res,
            comp_res_anno=comp_res_anno,
            df_for_comp=df_for_comp))

}

# drug_sen_df
# group_df
# drug_family_df
# g1,g2,g1_name,g2_name
# 
compare_drug_sen_two_drug <- function(comp_res,drug_family_df,g1,g2,g1_name,g2_name){
  for_family_test = comp_res%>%filter(p<0.05&group1==g1&group2==g2)

  for_family_test$which = NA
  for_family_test[for_family_test$mean_g1<for_family_test$mean_g2,]$which = g1_name
  for_family_test[for_family_test$mean_g1>for_family_test$mean_g2,]$which = g2_name

  for_family_test = merge(for_family_test,drug_family_df,by='inhibitor',all.x=TRUE)

  return(for_family_test)
}

if(FALSE){

  drug_family_df = read_tsv('~/my_project/data/VIZOME/nature_aml_drug_families.txt')
  drug_anno_df = read_tsv('~/my_project/data/VIZOME/nature_aml_drug_annotation.txt')
  drug_sen_df_viz = read_tsv('~/my_project/data/VIZOME/nature_aml_drug_sen.tsv')
  drug_sen_df_viz$lab_id = paste('X',str_replace(drug_sen_df_viz$lab_id,'-','.'),sep='')
  drug_sen_df_viz$sample = drug_sen_df_viz$lab_id

  compare_drug_sen_res = compare_drug_sen(drug_sen_df=drug_sen_df_viz, 
                                    group_df=subtyping_res$sample_group_df, 
                                    drug_anno_df=drug_anno_df, ic50auc='i50')

  compare_drug_sen_two_drug_res = compare_drug_sen_two_drug(comp_res=compare_drug_sen_res$comp_res,
                            drug_family_df=drug_family_df,
                            group1='c1',group2='c1',
                            group1_name='g1 carbohydrate sen',group2_name='g2 lipoid sen')

}


##############################################
# NetBID分析
##############################################


##############################################
# 转录因子分析
##############################################


##############################################
# 免疫组库分析 TRUST
##############################################







