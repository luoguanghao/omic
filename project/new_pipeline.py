
############################
############################
# Pipeline   #
############################
############################

############################
## analysis of GDSC+CCLE
############################
# raw tpm file from web
tpm_path = '/y/home/lgh/downloads/data/CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt'  ### <-
tpm_exp_mat = pd.read_csv(tpm_path, sep='\t')
%R library("AnnotationDbi")
%R library("org.Hs.eg.db")
ENSEMBL = np.array([i.split('.')[0] for i in tpm_exp_mat['gene_id']])
%R -i ENSEMBL
mapIds_ans = %R mapIds(org.Hs.eg.db,keys=ENSEMBL,column="SYMBOL",keytype="ENSEMBL",multiVals="first")

new_tpm_exp_mat = tpm_exp_mat.loc[mapIds_ans !='NA']
new_tpm_exp_mat.index = mapIds_ans[mapIds_ans !='NA']

new_tpm_exp_mat = new_tpm_exp_mat.drop(['gene_id','transcript_id'],axis=1)
new_tpm_exp_mat = new_tpm_exp_mat.loc[~new_tpm_exp_mat.index.duplicated(keep='first')]

new_tpm_exp_mat.to_csv('/y/Bondi/data/CCLE/expression/CCLE_RNAseq_tpm.tsv',sep='\t')

# good tpm file start from there
tpm_exp_mat = pd.read_csv('/y/Bondi/data/CCLE/expression/CCLE_RNAseq_tpm.tsv',sep='\t', index_col=0)
tpm_exp_mat.columns = [s.split('_')[0] for s in tpm_exp_mat.columns]
tpm_exp_mat = np.log(0.1+tpm_exp_mat)

## CCLE annotation
CCLE_annotation_file_path = '/y/Bondi/data/CCLE/sample_info/Cell_lines_annotations_20181226.txt'

CCLE_annotation_df = pd.read_csv(CCLE_annotation_file_path, sep='\t')
CCLE_annotation_df['CCLE_ID'] = [s.split('_')[0] for s in CCLE_annotation_df['CCLE_ID']] # 截短cell line名字

NAME_CCLEID_annotation_df = CCLE_annotation_df.loc[CCLE_annotation_df['type_refined'] == CCLE_annotation_df['type_refined']].set_index('Name')
CCLEID_NAME_annotation_df = CCLE_annotation_df.loc[CCLE_annotation_df['type_refined'] == CCLE_annotation_df['type_refined']].set_index('CCLE_ID')


# get  sp drug
drug = 'LFM-A13'

drug_sen_path = '/y/Bondi/data/GDSC/cell_line_sen/LFM-A13_DLBC.csv'

drug_sen_df = pd.read_csv(drug_sen_path).set_index('Cell line')
drug_sen_df = drug_sen_df.loc[set(NAME_CCLEID_annotation_df.index) & set(drug_sen_df.index)]
drug_sen_df.index = NAME_CCLEID_annotation_df.loc[drug_sen_df.index]['CCLE_ID']

drug_sen_series = drug_sen_df['IC50'].sort_values()
useful_sample = list(drug_sen_series.index)
sp_tpm_exp_mat = tpm_exp_mat[useful_sample]


# =================
# analysis sp drug
expr_mat = sp_tpm_exp_mat

sen_ls = list(drug_sen_series.iloc[:5].index)
unsen_ls = list(drug_sen_series.iloc[-5:].index)
condition_table=np.array([0]*5+[1]*5)

### 取top var
'''
expr_mat = expr_mat.loc[np.sum(expr_mat, axis=1) != 0]
top_var_gene = list(np.var(expr_mat, axis=1).sort_values().iloc[-5000:].index)
expr_mat = expr_mat.loc[top_var_gene]
'''
expr_mat = filter_by_var(expr_mat, which='gene', what='nozero')

p_res = plot_pheatmap_drug_sen(exp_mat=expr_mat.loc[np.var(expr_mat, axis=1).sort_values().iloc[-100:].index][sen_ls+unsen_ls], 
                               condition_table=np.array([0]*5+[1]*5))


# DESeq2
do_limma_res = do_limma(expr_mat=expr_mat.loc[:][sen_ls+unsen_ls], 
         condition_table=condition_table, 
         logFC_th=1, pvalue_th=0.05)

DEG_ls = do_limma_res['limma_output'].sort_values('P.Value').iloc[:100].index
p_res = plot_pheatmap_drug_sen(exp_mat=expr_mat.loc[DEG_ls][sen_ls+unsen_ls], 
                               condition_table=condition_table)

do_limma_res['limma_output'].to_csv('./%s_limma_df.tsv'%drug,sep='\t')

# select from deseq2 gene with use of RF
## candi_gene_ls = do_limma_res['limma_output'].sort_values('P.Value').iloc[:500].index
candi_gene_ls = np.array( do_limma_res['limma_output'].loc[(abs(do_limma_res['limma_output']['logFC'])>1) & (do_limma_res['limma_output']['P.Value']<0.05)].index )

do_random_forest_res = do_random_forest(X=expr_mat.loc[candi_gene_ls].T, 
                 y=np.log(np.array(drug_sen_series)), 
                 all_feature_ls=np.array(expr_mat.loc[candi_gene_ls].index), 
                 max_depth=2, 
                 random_state=random.randint(100,200), 
                 verbose=1)


rf_gene_ls = np.array(list(do_random_forest_res['feature_weight_df'].index))

p_res = plot_pheatmap_drug_sen(exp_mat=expr_mat.loc[rf_gene_ls][sen_ls+unsen_ls], 
                               condition_table=condition_table)

# EnRichment
enrich_result_df = do_enrichplot(rf_gene_ls)

do_deseq2_df = do_limma_res['limma_output']
deseq_feature_df = pd.DataFrame( {'ID':list(do_deseq2_df.index),
              'logFC':list(do_deseq2_df['log2FoldChange'])} )
prepare_for_GOcircle(feature_sorted_df=deseq_feature_df,
                    enrich_result_df=enrich_result_df, 
                    method='DESeq2_RF', 
                    save_path='./')







########################################################
## analysis of CH : DESeq2 , stepwise from DEG , EN , RF
########################################################

# import data
## drug sen
ch_drug_sen_path = '/y/Bondi/data/changhai_data/hl3-5/PDC体外用药数据汇总：包含AML和ALL.xlsx'
ch_drug_sen_df = load_changhai_drug_sen(ch_drug_sen_path)
drug_sen_series_1 = ch_drug_sen_df['PY34'].loc[ch_drug_sen_df['PY34'].isna()==False].sort_values()

drug_sen_df_2 = pd.read_csv('/y/Bondi/data/changhai_data/py34_2021-3-16-IC50新样本（这批样本可能没有测rnaseq，可能仅作参考）.txt',
                          sep='\t',
                          encoding='gb18030',
                          index_col=1,
                          header=1)
drug_sen_series_2 = drug_sen_df_2.loc['PY34'].iloc[3:].apply(clean_drug_value).sort_values()
drug_sen_series_2.index = [name[4:] for name in drug_sen_series_2.index]

drug_sen_series_all = pd.concat([drug_sen_series_1,drug_sen_series_2]).sort_values()
drug_sen_series_all = drug_sen_series_all[~drug_sen_series_all.index.duplicated(keep='first')]


## exp mat
count_path = '/y/Bondi/data/changhai_data/hl3-5/转录组原始数据（急性髓系白血病）：PDCallcount-onlyAML.csv'
count_expr_mat = pd.read_csv(count_path, sep=',', index_col=0)

count_path_2 = '/y/Bondi/data/changhai_data/hl3-5/ALL-10rnaseq.csv'
count_expr_mat2 = pd.read_csv(count_path_2, sep=',', index_col=0)

count_expr_mat_all = pd.concat([count_expr_mat,count_expr_mat2],axis=1)

# clean data
# 长海貌似不用clean

# get useful

useful_sample = list(set(drug_sen_series_all.index)&set(count_expr_mat_all.columns))
drug_sen_series = drug_sen_series_all.loc[useful_sample].sort_values()
count_expr_mat = count_expr_mat_all[useful_sample]

# =======================
# get sen unsen group

sen_ls = list(drug_sen_series.iloc[:10].index)
unsen_ls = list(drug_sen_series.iloc[-10:].index)
middle_ls = list(drug_sen_series.iloc[30:45].index)

# =========================
# filter by var and zeros

## 去掉全0的基因 取top var gene
count_expr_mat = count_expr_mat.loc[np.sum(count_expr_mat, axis=1) != 0]
top_var_gene = list(np.var(count_expr_mat, axis=1).sort_values().iloc[-5000:].index)
## 取top var
top100_var_gene = list(np.var(count_expr_mat, axis=1).sort_values().iloc[-100:].index)

# =========================================================
# Normalization and remove batch effect with use of DESeq2

condition_table = np.array([1]*int(count_expr_mat.shape[1]/2)+[0]*(count_expr_mat.shape[1]-int(count_expr_mat.shape[1]/2)))
do_deseq2_for_normalize = do_deseq2(expr_mat=count_expr_mat.loc[top_var_gene], 
                                    condition_table=condition_table, 
                                    colData=None)
## save
do_deseq2_for_normalize['expr_norm_mat'].to_csv('/y/Bondi/data/changhai_data/hl3-5/AML_ALL-10rnaseq_norm_deseq2.scv')
## load
exp_mat_deseq2_norm = pd.read_csv('/y/Bondi/data/changhai_data/hl3-5/AML_ALL-10rnaseq_norm_deseq2.scv')
## check with use of pheatmap
p_res = plot_pheatmap_drug_sen(exp_mat=filter_by_var(exp_mat_deseq2_norm[sen_ls+unsen_ls],what=1000), 
                               drug_sen_series=drug_sen_series,
                               condition_table=condition_table,
                               scale='column')

# ================================
# Feature Selection

# DESeq2

condition_table = np.array([0]*len(sen_ls)+[1]*len(unsen_ls))
res_do_deseq2 = do_deseq2(expr_mat=count_expr_mat.loc[top_var_gene][sen_ls+unsen_ls], 
                          condition_table=condition_table, colData=None)


# 逐步选择法选基因

condition_table = np.array([0]*len(sen_ls)+[1]*len(unsen_ls))
do_deseq2_df = pd.read_csv('/y/Bondi/jupyter/hl/3-23-PY34/result/py34_deseq2_df.tsv',sep='\t').set_index('gene_id') # get deseq2 result
exp_mat_deseq2_norm = pd.read_csv('/y/Bondi/data/changhai_data/hl3-5/AML_ALL-10rnaseq_norm_deseq2.scv').set_index('gene_id') # 

## 先看看前20个gene的分类能力
DEG_ls = do_deseq2_df.sort_values('padj').iloc[:20].index
p_res = plot_pheatmap_drug_sen(exp_mat=exp_mat_deseq2_norm.loc[DEG_ls][sen_ls+unsen_ls], 
                               condition_table=condition_table)

## 进行逐步选择 冗余
candi_gene_ls=np.array( do_deseq2_df.loc[(abs(do_deseq2_df['log2FoldChange'])>1) & (do_deseq2_df['pvalue']<0.05)].index ) # 取统计显著的
# DEG_ls = do_deseq2_df.sort_values('padj').iloc[:500].index # 取pvalue前500的
#len(DEG_gene_ls)
selected_feature_ls = step_wise_select_by_pheatmap(exp_mat=exp_mat_deseq2_norm.loc[candi_gene_ls][sen_ls+unsen_ls], 
                                                   condition_table=condition_table, 
                                                   candi_len=len(candi_gene_ls))

p_res = plot_pheatmap_drug_sen(exp_mat=exp_mat_deseq2_norm.loc[selected_feature_ls][sen_ls+unsen_ls], 
                               condition_table=condition_table)

## 进行逐步选择 冗余
def fry_selection():
    # DEG_gene_ls = do_deseq2_df.sort_values('padj').iloc[:500].index
    condition_table = [0]*len(sen_ls)+[1]*len(unsen_ls)
    tmp_feature_ls = [DEG_gene_ls[0]]
    mi = 0
    for g in DEG_gene_ls[1:500]:
        tmp_feature_ls.append(g)
        # mat_for_plot = res_do_deseq2['expr_norm_mat'].loc[tmp_feature_ls].T
        mat_for_plot = exp_mat_deseq2_norm[sen_ls+unsen_ls].loc[tmp_feature_ls].T
        %R library(pheatmap)
        %R -i mat_for_plot
        %R p_res = pheatmap(mat_for_plot,scale='column',cutree_row = 2,clustering_method = "ward.D",clustering_distance_rows = "correlation",silent=TRUE)
        # order = %R p_res$tree_row$order
        # labels = %R p_res$tree_row$label
        row_cluster = %R cutree(p_res$tree_row,k=2)
        tmp_mi = metrics.adjusted_mutual_info_score(condition_table, row_cluster)
        #print('mi now is %s'%tmp_mi)
        if tmp_mi <= mi:
            #print('drop')
            tmp_feature_ls = tmp_feature_ls[:-1]
        else:
            mi = tmp_mi
    print(mi)
    fry_feature_ls = tmp_feature_ls


# Elastic Net选择
candi_gene_ls = do_deseq2_df.sort_values('padj').iloc[:500].index # 候选gene list的选择，可以是DEG top gene，也可以是别的基因
res_do_en = do_elastic_net(X=np.log(exp_mat_deseq2_norm.loc[candi_gene_ls].T), 
                            y=np.log(np.array(drug_sen_series)), 
                            gene_ls=np.array(exp_mat_deseq2_norm.loc[candi_gene_ls].index), 
                            alphas=None, 
                            cv=10, 
                            max_iter=5000, 
                            random_state=random.randint(100,200), 
                            selection='random', 
                            verbose=1)
en_gene_ls = res_do_en['feature_ls']

# random forest 选择
candi_gene_ls = do_deseq2_df.sort_values('padj').iloc[:500].index
do_random_forest_res = do_random_forest(X=np.log(exp_mat_deseq2_norm.loc[candi_gene_ls].T), 
                 y=np.log(np.array(drug_sen_series)), 
                 all_feature_ls=np.array(exp_mat_deseq2_norm.loc[candi_gene_ls].index), 
                 max_depth=2, 
                 random_state=random.randint(100,200), 
                 verbose=1)
rf_gene_ls = do_random_forest_res['feature_ls']

p_res = plot_pheatmap_drug_sen(exp_mat=exp_mat_deseq2_norm.loc[rf_gene_ls][sen_ls+unsen_ls], 
                               condition_table=condition_table)

# WGCNA
condition_table = np.array([0]*len(sen_ls)+[1]*len(unsen_ls))
exp_mat_deseq2_norm = pd.read_csv('/y/Bondi/data/changhai_data/hl3-5/AML_ALL-10rnaseq_norm_deseq2.scv')

res_WGCNA_analysis = WGCNA_analysis(expr_mat=exp_mat_deseq2_norm[sen_ls+unsen_ls], 
                                    drug_sen_series=drug_sen_series[sen_ls+unsen_ls], 
                                    n_top_mad=2000, 
                                    plot_flag=True)
res_WGCNA_analysis['moduleEigengenes_drug_cor_df']
wgcna_gene = list(res_WGCNA_analysis['cluster_df'].loc['green']['gene'])
p_res = plot_pheatmap_drug_sen(exp_mat=exp_mat_deseq2_norm.loc[wgcna_gene][sen_ls+unsen_ls], 
                               condition_table=condition_table)

# NetBID2



# ================================
# POST Feature Selection

# 各种 feature 画热图
deq2_DEG = np.array( 
    do_deseq2_df.loc[(abs(do_deseq2_df['log2FoldChange'])>1) & (do_deseq2_df['pvalue']<0.05)].index )

p_res = plot_pheatmap_drug_sen(exp_mat=exp_mat_deseq2_norm.loc[deq2_DEG][sen_ls+unsen_ls], 
                               condition_table=condition_table)


# 和tf或者sig做交集

%R library(NetBID2)
%R db.preload(use_level='gene',use_spe='human',update=FALSE)
tf_ls = %R tf_sigs$tf$external_gene_name
sig_ls = %R tf_sigs$sig$external_gene_name
overlap_feature_ls = list(set(deq2_DEG)&set(tf_ls))
p_res = plot_pheatmap_drug_sen(exp_mat=exp_mat_deseq2_norm.loc[set(deq2_DEG)&set(tf_ls)][sen_ls+unsen_ls], 
                               condition_table=condition_table)

# EnRich & GOcircle

enrich_result_df = do_enrichplot(rf_feature_array)

deseq_feature_df = pd.DataFrame( {'ID':list(do_deseq2_df.index),
              'logFC':list(do_deseq2_df['log2FoldChange'])} )
prepare_for_GOcircle(feature_sorted_df=deseq_feature_df,
                    enrich_result_df=enrich_result_df, 
                    method='DESeq2', 
                    save_path='./')




###################################################
## analysis of VIZOME : limma , WGCNA , cor one by one
## constrict FLT-ITD sample
###################################################
drug_ls = ['Midostaurin', 'Quizartinib (AC220)', 'Crenolanib', 'Gilteritinib (ASP-2215)']

# import data
drug_sen_path = '/y/Bondi/data/VIZOME/nature_aml_drug_sen.tsv'
drug_sen_df = pd.read_csv(drug_sen_path, sep='\t', index_col=0)

rpkm_path = '/y/Bondi/data/VIZOME/nature_aml_log2_fpkm.txt'
rpkm_expr_mat = pd.read_csv(rpkm_path, sep='\t', index_col=0)

Variants_path =  '/y/Bondi/data/VIZOME/LAML_nature_PDC.maf'
Variants_df = pd.read_csv(Variants_path, sep='\t')
## 获取 FLT3-ITD 的病人
FLT3ITD_Variants_df = Variants_df.loc[Variants_df['Hugo_Symbol']=='FLT3'].loc[Variants_df['ITDorNOT']=='ITD']
FLT3ITD_sample_ls = list(set(FLT3ITD_Variants_df['Tumor_Sample_Barcode']))

## get FLT-ITD sample
rpkm_expr_mat = rpkm_expr_mat[~rpkm_expr_mat.index.duplicated(keep='first')]
FLT3ITD_useful_sample_ls = list(set(rpkm_expr_mat.columns) & set(FLT3ITD_sample_ls))


# filter gene
### 取top var
'''
expr_mat = rpkm_expr_mat
expr_mat = expr_mat.loc[np.sum(expr_mat, axis=1) != 0]
top_var_gene = list(np.var(expr_mat, axis=1).sort_values().iloc[-5000:].index)
expr_mat = expr_mat.loc[top_var_gene]
rpkm_expr_mat = expr_mat
'''
rpkm_expr_mat = filter_by_var(rpkm_expr_mat, which='gene', what='nozero')

# DO
do_diff_corr_wgcna_res = do_diff_corr_wgcna(drug_ls[drug_id], 
                   all_drug_sen_df=all_drug_sen_df, 
                   rpkm_expr_mat=rpkm_expr_mat[FLT3ITD_useful_sample_ls], 
                   get_sen_unsen_p=0.2 , 
                   get_sen_unsen_c = None)

## Post DO
output_path = '#############'
## WGCNA 看看哪个module好 选哪个
do_diff_corr_wgcna_res['WGCNA_analysis_res']['moduleEigengenes_drug_cor_df'].sort_values('Adjusted R-squared')
wgcna_ls_0 = do_diff_corr_wgcna_res['WGCNA_analysis_res']['cluster_df'].loc['green']['gene']
wgcna_ls_0.to_csv('%s/do_WGCNA_gene_ls_%s.tsv'%(output_path, drug_ls[drug_id]), sep='\t')

## limma
do_diff_corr_wgcna_res['do_limma_res']['limma_output'].sort_values('adj.P.Val') # look look
do_diff_corr_wgcna_res['do_limma_res']['limma_output'].sort_values('adj.P.Val').to_csv(
                                            '%s/do_limma_res_%s.tsv'%(output_path, drug_ls[drug_id]), sep='\t')

## cor
do_diff_corr_wgcna_res['cor_res_df'].sort_values('pvalue') # look look
do_diff_corr_wgcna_res['cor_res_df'].sort_values('pvalue').to_csv(
                                            '%s/do_cor_res_%s.tsv'%(output_path, drug_ls[drug_id]), sep='\t')



#########################################
# TCGA & GTEx
#########################################

TCGA_exp_mat = pd.read_csv('/y/Bondi/data/TCGA/LAML/TCGA-LAML.htseq_counts.tsv.gz',sep='\t').set_index('Ensembl_ID')
id_symbol_Series = pd.read_csv('/y/Bondi/data/TCGA/LAML/gencode.v22.annotation.gene.probeMap',sep='\t').set_index('id')['gene']
TCGA_exp_mat = TCGA_data_clean(TCGA_exp_mat,id_symbol_Series)['exp_mat']

gct_path = '/y/Bondi/data/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz'
GTEx_exp_mat = load_GTEx_gct_file(gct_path, dup_kep='first')

Anno_df = pd.read_csv('/y/Bondi/data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.tsv',sep='\t').set_index('SAMPID')

blood_sample_gtex = set(Anno_df.loc[Anno_df['SMTS']=='Blood'].index)&set(GTEx_exp_mat.columns)

MVA_pathway_genes = list('ACAT1 ACAT2 HMGCS1 HMGCR MVK PMVK MVD IDI1 IDI2'.split(' '))

# combine TCGA GTEx
TCGA_GTEx_exp_mat = pd.concat([TCGA_exp_mat.loc[set(TCGA_exp_mat.index)&set(GTEx_exp_mat.index)],
           GTEx_exp_mat[blood_sample_gtex].loc[set(TCGA_exp_mat.index)&set(GTEx_exp_mat.index)]],axis=1).astype(int)

# TCGA_GTEx_exp_mat.to_csv('./TCGA_GTEx_exp_mat_count.tsv',sep='\t')

tcga_sam_n = 151
gtex_sam_n = 110
condition_table = np.array([0]*151+[1]*110)
# Normalization with use of deseq2
do_deseq2_res = do_deseq2(expr_mat=TCGA_GTEx_exp_mat, condition_table=condition_table, colData=None)

MVA_pathway_mean_exp_series = np.mean(do_deseq2_res['expr_norm_mat'].loc[MVA_pathway_genes])

expr_norm_mat = do_deseq2_res['expr_norm_mat']

df_for_plot = pd.DataFrame({
        'MVA_pathway_mean_exp':list(MVA_pathway_mean_exp_series),
        'level':['Tumor']*151+['Normal']*110 })
ggplot2_boxplot(df=df_for_plot, x_level='level', y_value='MVA_pathway_mean_exp', stat_method=np.nan)



############################
# barplot for mut & drug sen
############################

drug_sen_path = '/y/Bondi/data/VIZOME/nature_aml_drug_sen.tsv'
Variants_path = '/y/Bondi/data/VIZOME/LAML_nature_PDC.maf'
drug_sen_df = pd.read_csv(drug_sen_path, sep='\t', index_col=0)
Variants_df = pd.read_csv(Variants_path, sep='\t')

## 获取 FLT3-ITD 的病人
FLT3ITD_Variants_df = Variants_df.loc[Variants_df['Hugo_Symbol']=='FLT3'].loc[Variants_df['ITDorNOT']=='ITD']
FLT3ITD_sample_ls = list(set(FLT3ITD_Variants_df['Tumor_Sample_Barcode']))
FLT3ITD_sample_Variants_df = Variants_df.loc[[i in FLT3ITD_sample_ls for i in Variants_df['Tumor_Sample_Barcode']]]

from collections import Counter
def get_common_mut(Variants_df, n=10):
    sample_count = len(set(Variants_df['Tumor_Sample_Barcode']))
    
    mut_count_df = pd.DataFrame(
        Counter(
            Variants_df['Hugo_Symbol']).items(),columns=['gene','count']).set_index('gene').sort_values('count',ascending=False)    
    
    mut_count_df['frequency'] = mut_count_df['count']/sample_count
    #most_common_res = Counter(Variants_df['Hugo_Symbol']).most_common(10)
    #top_genes = [g[0] for g in most_common_res]
    top_gene_ls = list(mut_count_df.index[:n])
    return {'mut_count_df':mut_count_df,
           'top_gene_ls':top_gene_ls}

top_genes = get_common_mut(FLT3ITD_sample_Variants_df)['top_gene_ls']

drug = 'Gilteritinib (ASP-2215)'

drug_sen_series = drug_sen_df.loc[drug].set_index(['lab_id'])['auc']
# Variants_df
specific_sample_ls = FLT3ITD_sample_ls
plot_gene_ls = top_genes
text = '%s: with FLTIDT'%drug

bar_plot_mut_drug_sen(drug_sen_series,specific_sample_ls,plot_gene_ls,text)



















































































































































































































###########
# 取topvar，除去0的gene
###########

# 去掉全0的基因 取top var gene
count_expr_mat = count_expr_mat.loc[np.sum(count_expr_mat, axis=1) != 0]
top_var_gene = list(np.var(count_expr_mat, axis=1).sort_values().iloc[-5000:].index)

# 可选取top100 gene
# top100_var_gene = list(np.var(count_expr_mat, axis=1).sort_values().iloc[-100:].index)

###############
# 拿DESeq2对数据进行normalization
# 也可以直接获取已经normalization的数据
###############
condition_table = np.array([1]*int(count_expr_mat.shape[1]/2)+[0]*(count_expr_mat.shape[1]-int(count_expr_mat.shape[1]/2)))
do_deseq2_for_normalize = do_deseq2(expr_mat=count_expr_mat.loc[top_var_gene], 
                                    condition_table=condition_table, 
                                    colData=None)

# do_deseq2_for_normalize['expr_norm_mat'].to_csv('/y/Bondi/data/changhai_data/hl3-5/AML_ALL-10rnaseq_norm_deseq2.scv')

exp_mat_deseq2_norm = pd.read_csv('/y/Bondi/data/changhai_data/hl3-5/AML_ALL-10rnaseq_norm_deseq2.scv')


# check top 100 var 的分类效果
mat_for_plot = do_deseq2_for_normalize['expr_norm_mat'].loc[top100_var_gene][sen_ls+unsen_ls].T
%R library(pheatmap)
%R -i mat_for_plot
%R p_res = pheatmap(mat_for_plot,scale='column')
order = %R p_res$tree_row$order
labels = %R p_res$tree_row$label
plt.scatter(
    np.arange(len(labels)),
    drug_sen_series.loc[mat_for_plot.index].loc[labels[order-1]],
            )


##########
# deseq2
##########

condition_table = np.array([0]*len(sen_ls)+[1]*len(unsen_ls))
res_do_deseq2 = do_deseq2(expr_mat=count_expr_mat.loc[top_var_gene][sen_ls+unsen_ls], 
                          condition_table=condition_table, colData=None)

DEG_ls = res_do_deseq2['DESeq_res'].sort_values('padj').iloc[:20].index
mat_for_plot = res_do_deseq2['expr_norm_mat'].loc[DEG_ls].T
%R library(pheatmap)
%R -i mat_for_plot
%R p_res = pheatmap(mat_for_plot,scale='column')
order = %R p_res$tree_row$order
labels = %R p_res$tree_row$label
plt.scatter(
    np.arange(len(labels)),
    drug_sen_series.loc[mat_for_plot.index].loc[labels[order-1]],
            )


## DESeq2 逐步回归
do_deseq2_df = pd.read_csv('/y/Bondi/jupyter/hl/3-23/result/py34_deseq2_df.tsv',sep='\t').set_index('gene_id')
exp_mat_deseq2_norm = pd.read_csv('/y/Bondi/data/changhai_data/hl3-5/AML_ALL-10rnaseq_norm_deseq2.scv').set_index('gene_id')

DEG_gene_ls = do_deseq2_df.sort_values('padj').iloc[:500].index

condition_table = [0]*len(sen_ls)+[1]*len(unsen_ls)

tmp_feature_ls = [DEG_gene_ls[0]]
mi = 0
for g in DEG_gene_ls[1:500]:
    tmp_feature_ls.append(g)

    # mat_for_plot = res_do_deseq2['expr_norm_mat'].loc[tmp_feature_ls].T
    mat_for_plot = exp_mat_deseq2_norm[sen_ls+unsen_ls].loc[tmp_feature_ls].T
    %R library(pheatmap)
    %R -i mat_for_plot
    %R p_res = pheatmap(mat_for_plot,scale='column',cutree_row = 2,clustering_method = "ward.D",clustering_distance_rows = "correlation",silent=TRUE)
    # order = %R p_res$tree_row$order
    # labels = %R p_res$tree_row$label
    row_cluster = %R cutree(p_res$tree_row,k=2)
    tmp_mi = metrics.adjusted_mutual_info_score(condition_table, row_cluster)
    #print('mi now is %s'%tmp_mi)
    if tmp_mi < mi:
        #print('drop')
        tmp_feature_ls = tmp_feature_ls[:-1]
    else:
        mi = tmp_mi

print(mi)



def step_wise_select_by_pheatmap(exp_mat, condition_table):
    '''
    exp_mat: gene*sample
    condition_table: as long as exp_mat.columns
    
    exp_mat = exp_mat_deseq2_norm[sen_ls+unsen_ls].loc[DEG_gene_ls]
    condition_table = [0]*len(sen_ls)+[1]*len(unsen_ls)
    '''
    # DEG_gene_ls = do_deseq2_df.sort_values('padj').iloc[:500].index
    # condition_table = [0]*len(sen_ls)+[1]*len(unsen_ls)
    candi_feature_ls = list(exp_mat.columns)

    tmp_feature_ls = [candi_feature_ls[0]]
    mi = 0
    for g in candi_feature_ls[1:500]:
        tmp_feature_ls.append(g)

        # mat_for_plot = res_do_deseq2['expr_norm_mat'].loc[tmp_feature_ls].T
        mat_for_plot = exp_mat.loc[tmp_feature_ls].T
        %R library(pheatmap)
        %R -i mat_for_plot
        %R p_res = pheatmap(mat_for_plot,scale='column',cutree_row = 2,clustering_method = "ward.D",clustering_distance_rows = "correlation",silent=TRUE)
        # order = %R p_res$tree_row$order
        # labels = %R p_res$tree_row$label
        row_cluster = %R cutree(p_res$tree_row,k=2)
        tmp_mi = metrics.adjusted_mutual_info_score(condition_table, row_cluster)
        #print('mi now is %s'%tmp_mi)
        if tmp_mi < mi:
            #print('drop')
            tmp_feature_ls = tmp_feature_ls[:-1]
        else:
            mi = tmp_mi

    print(mi)

    selected_feature_ls = tmp_feature_ls
    return selected_feature_ls



def plot_pheatmap_drug_sen(exp_mat, condition_table):
    '''
    exp_mat : gene*sample
    condition_table : the same as len of sample

    exp_mat = exp_mat_deseq2_norm[sen_ls+unsen_ls].loc[tmp_feature_ls]
    '''
    mat_for_plot = exp_mat.T
    %R library(pheatmap)
    %R -i mat_for_plot
    %R p_res = pheatmap(mat_for_plot,scale='column',cutree_row = 2,clustering_method = "ward.D",clustering_distance_rows = "correlation")
    order = %R p_res$tree_row$order
    labels = %R p_res$tree_row$label
    plt.scatter(
        np.arange(len(labels)),
        drug_sen_series.loc[mat_for_plot.index].loc[labels[order-1]],
                )

    row_cluster = %R cutree(p_res$tree_row,k=2)
    mi = metrics.adjusted_mutual_info_score(condition_table, row_cluster)

    #return {'group_df':group_df,
    #        'mi':mi
    #        }

    # DEG_ls = res_do_deseq2['DESeq_res'].sort_values('padj').iloc[:100].index
    mat_for_plot = exp_mat_deseq2_norm[sen_ls+unsen_ls].loc[tmp_feature_ls].T
    %R library(pheatmap)
    %R -i mat_for_plot
    %R p_res = pheatmap(mat_for_plot,scale='column',cutree_row = 2,clustering_method = "ward.D",clustering_distance_rows = "correlation")
    order = %R p_res$tree_row$order
    labels = %R p_res$tree_row$label
    plt.scatter(
        np.arange(len(labels)),
        drug_sen_series.loc[mat_for_plot.index].loc[labels[order-1]],
                )


############
# Elastic Net
############
DEG_gene_ls = do_deseq2_df.sort_values('padj').iloc[:500].index
res_do_en = do_elastic_net(X=np.log(exp_mat_deseq2_norm.loc[DEG_gene_ls].T), 
                                  y=np.log(np.array(drug_sen_series)), 
                                  gene_ls=np.array(exp_mat_deseq2_norm.loc[DEG_gene_ls].index), 
                                  alphas=None, 
                                  cv=10, 
                                  max_iter=5000, 
                                  random_state=random.randint(100,200), 
                                  selection='random', 
                                  verbose=1)


en_gene_ls = list(res_do_en['feature_weight_df'].index)

mat_for_plot = np.log(exp_mat_deseq2_norm.loc[en_gene_ls][sen_ls+unsen_ls].T)
%R library(pheatmap)
%R -i mat_for_plot
%R p_res = pheatmap(mat_for_plot, scale='column')

order = %R p_res$tree_row$order
labels = %R p_res$tree_row$label
plt.scatter(
    np.arange(len(labels)),
    drug_sen_series.loc[mat_for_plot.index].loc[labels[order-1]],
            )


##################
# Random Forest
##################
DEG_gene_ls = do_deseq2_df.sort_values('padj').iloc[:500].index
do_random_forest_res = do_random_forest(X=np.log(exp_mat_deseq2_norm.loc[DEG_gene_ls].T), 
                 y=np.log(np.array(drug_sen_series)), 
                 all_feature_ls=np.array(exp_mat_deseq2_norm.loc[DEG_gene_ls].index), 
                 max_depth=2, 
                 random_state=random.randint(100,200), 
                 verbose=1)


rf_gene_ls = list(do_random_forest_res['feature_weight_df'].index)

mat_for_plot = np.log(exp_mat_deseq2_norm.loc[rf_gene_ls][sen_ls+unsen_ls].T)
%R library(pheatmap)
%R -i mat_for_plot
%R p_res = pheatmap(mat_for_plot, scale='column')

order = %R p_res$tree_row$order
labels = %R p_res$tree_row$label
plt.scatter(
    np.arange(len(labels)),
    drug_sen_series.loc[mat_for_plot.index].loc[labels[order-1]],
            )








############
############
# 八卦图 GORich
############
############

deq2_DEG = np.array( do_deseq2_df.loc[(abs(do_deseq2_df['log2FoldChange'])>1) & (do_deseq2_df['pvalue']<0.05)].index )

deq2_DEG_df = pd.DataFrame({'ID':do_deseq2_df.loc[deq2_DEG]['log2FoldChange'].index,
              'logFC':list(do_deseq2_df.loc[deq2_DEG]['log2FoldChange'])})
deq2_DEG_df.to_csv('/y/Bondi/jupyter/hl/3-23/deq2_DEG_df_for_GOPlot.csv')


%R library("clusterProfiler")
%R library("org.Hs.eg.db")
%R library("enrichplot")
%R library("ggplot2")


'''
- 文件名
- Diff 基因列表
'''
DEG_gene_array = deq2_DEG # maybe other

%R -i DEG_gene_array
# %R DEG_gene_array <- deq2_DEG
%R df2 <- bitr(DEG_gene_array,fromType = "SYMBOL",toType = c("ENTREZID", "GENENAME"),OrgDb = org.Hs.eg.db)
DEG_gene_eid = %R df2[[2]]

%R -i DEG_gene_eid
%R kk=enrichGO(gene=DEG_gene_eid,OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, ont="ALL", readable=T)

enrich_result_df = %R as.data.frame(kk)
enrich_result_df = enrich_result_df[['ONTOLOGY','ID','Description','geneID','p.adjust']]
enrich_result_df.columns = ['Category','ID','Term','Genes','adj_pval']

enrich_result_df['Genes'] = enrich_result_df['Genes'].map(lambda x: x.replace('/',','))

enrich_result_df.index = list(range(1,len(enrich_result_df)+1))   ###

enrich_result_df.to_csv('./enrich_result_for_GOcircle_deq2_DEG_df.tsv',sep='\t') # maybe other





























