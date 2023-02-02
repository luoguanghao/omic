

from sklearn import linear_model
from sklearn.decomposition import PCA


def calculate_avg_score(cal_score_expr_mat):
    score = cal_score_expr_mat.mean(axis=0)
    # find the cutoff


def calculate_ssgsea_score(expr_mat, feature_gene_ls):
    gsva_gene_df = {'gene': feature_gene_ls,
                    'set_name': ['my_feature']*len(feature_gene_ls)}
    gsva_matrix = do_GSVA(expr_mat=expr_mat,
                          gene_df=gsva_gene_df)


def calculate_pca_score():
    pca = PCA(n_components=2)
    pca.fit(X)
    pca_score = pca.transform(X)


def calculate_linear_model_score(cal_score_expr_mat, drug_sen_Series):
    '''
    cal_score_expr_mat: genes × samples
    待续
    '''
    X = np.array(cal_score_expr_mat.T)
    reg = linear_model.RidgeCV(alphas=np.logspace(-6, 6, 13))
    reg.fit(X, [0, .1, 1])


def do_cal_score(expr_mat):
    # import data
    expr_mat_path = '/y/home/lgh/lgh/work/drug_sensitivity/data/all_breast_cell_com_geneexpression_norm.csv'
    drug_response_path = '/y/home/lgh/lgh/work/drug_sensitivity/data/ic50_formatted.csv'

    import_data_res = import_data(genomic_alter_path=None,
                                  expr_mat_path=expr_mat_path, drug_response_path=drug_response_path)
    expr_mat = import_data_res['expr_mat']
    drug_response_df = import_data_res['drug_response_df']

    depmap_info_info = '/y/Bondi/data/depmap/info/sample_info.csv'
    depmap_info_df = pd.read_csv(depmap_info_info)

    rpkm_CCLE_path = '/y/Bondi/data/CCLE/expression/CCLE_RNAseq_genes_rpkm_20180929.gct'
    CCLE_annotation_file_path = '/y/Bondi/data/CCLE/sample_info/Cell_lines_annotations_20181226.txt'
    rpkm_CCLE_df = utilities.file_tools.read_txt(
        rpkm_CCLE_path, columns_name_line=2)
    CCLE_annotation_df = pd.read_csv(CCLE_annotation_file_path, sep='\t')

    rpkm_CCLE_expr_mat = rpkm_CCLE_df.drop(
        labels='Name', axis=1).set_index('Description').apply(pd.to_numeric)
    rpkm_CCLE_expr_mat = clean_data(rpkm_CCLE_expr_mat)

    # 把CCLE expression的cell line name 转变为普通的line name
    # some cell line should be drop out
    rpkm_CCLE_expr_mat = rpkm_CCLE_expr_mat.drop(
        columns=['KE97_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 'NCIH684_LIVER', 'NHAHTDD_CENTRAL_NERVOUS_SYSTEM', 'SF767_CERVIX'])
    CCLE_cell_name_ls = list(rpkm_CCLE_expr_mat.columns)
    ccleName2cellLineName = depmap_info_df[[
        'CCLE_Name', 'cell_line_name']].set_index(['CCLE_Name'])
    rpkm_CCLE_expr_mat.columns = list(
        ccleName2cellLineName.loc[CCLE_cell_name_ls[:]]['cell_line_name'])

    #####
    # 导入GDSC数据

    ####
    AZD7762_GDSC_brca_df = pd.read_csv(
        '/y/Bondi/data/GDSC/AZD7762_GDSC1.tab', sep='\t').set_index(['Cell line'])
    AZD7762_GDSC_brca_cell_set = set(AZD7762_GDSC_brca_df.index)

    # ssgsea score
    feature_gene_ls = Kai_Feature
    gsva_gene_df = pd.DataFrame({
        'gene': feature_gene_ls,
        'set_name': ['my_feature']*len(feature_gene_ls)})

    gsva_gene_df = pd.DataFrame({'gene': feature_gene_ls, 'set_name': [
                                'my_feature']*len(feature_gene_ls)})
    gsva_matrix = do_GSVA(expr_mat=rpkm_CCLE_expr_mat[AZD7762_brca_CCLE_GDSC_cell_set],
                          gene_df=gsva_gene_df)

    y = np.array(
        AZD7762_GDSC_brca_df.loc[AZD7762_brca_CCLE_GDSC_cell_set]['IC50'])
    x = np.array(gsva_matrix.T['my_feature'])

    %R library(maxstat)
    %R - i x
    %R - i y
    %R mydata < - data.frame(cbind(x, y))
    %R mod < - maxstat.test(y ~ x, data=mydata, smethod="Wilcoxon", pmethod="HL", minprop=0.25, maxprop=0.75, alpha=0.05)
    %R print(mod)

    cutoff = %R mod$estimate
    cutoff = cutoff[0]

    group_low_score = []
    group_high_score = []

    for i in range(len(x)):
        if x[i] > cutpoint:
            group_high_score.append(y[i])
        else:
            group_low_score.append(y[i])

    plt.boxplot([group_low_score, group_high_score])
