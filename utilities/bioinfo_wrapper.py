import gseapy

########
# gsea #
########

def load_gmt_file(gmt_path, set_name_col=0, gene_name_col=1):
    # gmt_path = '/y/Bondi/data/pathway/reactome.gmt'
    file = open(gmt_path).read().strip().split('\n')
    pathway_gene_dict = {}
    for line in file:
        list_ls = line.strip().split('\t')
        pathway_gene_dict[list_ls[set_name_col]] = list_ls[gene_name_col:]

    return pathway_gene_dict

def gmt_2_do_GSVA_gene_set():
    pass


import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
def do_GSVA(expr_mat, gene_df):
    '''
    expr_mat: genes × samples Normalization!!!!
    gene_df: df:
        gene   set_name
        --------------
        gene1  A
        gene2  A
        gene3  B
    '''
    pandas2ri.activate() ## 对于rpy2转pandas很有用

    with localconverter(ro.default_converter + pandas2ri.converter):
        r_expr_mat = ro.conversion.py2rpy(expr_mat)
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_gene_df = ro.conversion.py2rpy(gene_df)

    rcode = """
    library(GSVA)
    gene_df <- %s
    gene_df <- as.data.frame(gene_df)
    expr_mat <- %s

    list<- split(as.matrix(gene_df)[,1], gene_df[,2])
    gsva_matrix<- gsva(as.matrix(expr_mat), list,method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE, parallel.sz=1)
    gsva_matrix <- as.data.frame(gsva_matrix)
    """%(r_gene_df.r_repr(), r_expr_mat.r_repr())

    ro.r(rcode)

    gsva_matrix = ro.r("gsva_matrix")
    return gsva_matrix


    '''
    %R -i expr_mat
    %R -i set_info_dfmat: gene × sample 要是预处理过的mat:除去重复行，滤掉无用genes
    gene_set: df:
            gene   set_name
            --------------
            gene1  A
            gene2  A
    %R list<- split(as.matrix(mmc_df)[,1], mmc_df[,2])
    =
    %R gsva_matrix<- gsva(as.matrix(rpkm_expr_mat), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

    '''


def ssGSEA_score(feature_gene_ls, expr_mat):
    '''
    feature_gene_ls: a list of genes
    expr_mat: maybe rpkm_expr_mat, no log2
    '''

    gsva_gene_df = pd.DataFrame({
        'gene': feature_gene_ls,
        'set_name': ['my_feature']*len(feature_gene_ls)})

    gsva_gene_df = pd.DataFrame({'gene': feature_gene_ls, 'set_name': [
                                'my_feature']*len(feature_gene_ls)})
    ssGSEA_score_series = do_GSVA(expr_mat=expr_mat,
                        gene_df=gsva_gene_df)

    return ssGSEA_score_series


def get_pathway_expression():
    pass


















































