



def clean_data(expr_mat, dup='avg'):
    '''
    Remove duplicated row according to their index and columns content  
    Remove the abnormal sample or gene  
    input:  
    - annotation dataframe and expression matrix
    - dup = avg/plus/first/discard
    return: a cleaned expr_mat  
    '''
    # deal with duplicated
    no_duplicated_expr_mat = None
    duplicated_expr_mat = None

    # gene_list = df_annotation.loc[expr_mat.index,'gene']
    sample_list = expr_mat.columns

    dup_row = expr_mat.index.duplicated(keep=False)
    no_duplicated_expr_mat = expr_mat[~dup_row]
    duplicated_expr_mat = expr_mat[dup_row]
    
    # 对重复基因提供四种处理策略:平均值,相加,第一值,删掉
    if dup == 'avg':
        ## avg:
        duplicated_gene_set_ls = list(set(duplicated_expr_mat.index))
        avg_duplicated_expr_mat = np.zeros([len(duplicated_gene_set_ls), len(sample_list)])
        for i_g in range(len(duplicated_gene_set_ls)):
            for i_s in range(len(sample_list)):
                avg_duplicated_expr_mat[i_g][i_s] = np.mean(
                    duplicated_expr_mat.loc[duplicated_gene_set_ls[i_g]][sample_list[i_s]]
                    )
        
        avg_duplicated_expr_mat = pd.DataFrame(
            avg_duplicated_expr_mat, 
            index=duplicated_gene_set_ls, 
            columns=sample_list)
    
        unduplicated_expr_mat = pd.concat(
            [no_duplicated_expr_mat, avg_duplicated_expr_mat], axis=0)

    elif dup == 'plus':    
        ## plus
        pass
    
    elif dup == 'plus':  
        ## first
        pass

    elif dup == 'plus':  
        ## discard
        pass

    return expr_mat[~expr_mat.index.duplicated(keep='first')]



def load_TCGA_data(expr_mat_path, annotation_path, survival_path):
    # load data
    ## expr
    # expr_mat_path = '/y/Bondi/data/TCGA/BRCA/TCGA-BRCA.htseq_fpkm.tsv.gz'
    expr_mat = pd.read_csv(expr_mat_path, sep='\t').set_index('Ensembl_ID')
  
    ## annotation
    # annotation_path = '/y/Bondi/data/TCGA/BRCA/gencode.v22.annotation.gene.probeMap'
    annotation_df = pd.read_csv(annotation_path, sep='\t').set_index('id')
  
    ## survival
    # survival_path = '/y/Bondi/data/TCGA/BRCA/TCGA-BRCA.survival.tsv.gz'
    survival_df = pd.read_csv(survival_path, sep='\t').set_index('sample')
  
    ## other



    # transform
    gene_ls = list(annotation_df.loc[expr_mat.index]['gene'])
    gene_ls = np.array([i.replace('-','_') for i in gene_ls]) # raw gene ls
    expr_mat.index = gene_ls

    # get useful sample
    sample_ls = np.array( list(set(expr_mat.columns)&set(survival_df.index)) )
    expr_mat = expr_mat[sample_ls]
    survival_df = survival_df.loc[sample_ls]
    
    sample_ls = np.array([i.replace('-','_') for i in sample_ls])
    expr_mat.columns = sample_ls
    survival_df.index = sample_ls

    # discart 多余 genes
    expr_mat = clean_data(expr_mat)

    gene_ls = np.array(expr_mat.index)


    return {
        'sample_ls':sample_ls,
        'gene_ls':gene_ls,
        'expr_mat':expr_mat,
        'survival_df':survival_df
    }







