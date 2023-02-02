import pandas as pd
import numpy as np




##########################
# read count 2 rpkm
##########################

def merge(intervals):
    intervals.sort(key=lambda x: x[0])

    merged = []
    for interval in intervals:
        # 如果列表为空，或者当前区间与上一区间不重合，直接添加
        if not merged or merged[-1][1] < interval[0]:
            merged.append(interval)
        else:
            # 否则的话，我们就可以与上一区间进行合并
            merged[-1][1] = max(merged[-1][1], interval[1])

    return merged


def readCount2Rpkm(gtf_path, expr_mat):
    '''
    import:
        gtf
        expr_mat: gene*sample
    output:
        normalization expr_mat : log2 rpkm
    '''
    # load gtf file
    f1 = open(gtf_path).read().strip().split('\n')
    gtf_dict = {'type':[], 'start':[],'end':[],'gene_name':[]}
    cnt= 0 
    for line in f1:
        if line[0]=='#':
            continue
        cnt += 1
        #if cnt > 10:
        #   break
        line_ls = line.split('\t')
        #print(line_ls[2],line_ls[3],line_ls[4], line.split('\t')[-1].strip().split(';')[0][9:-1])
        if line_ls[2]== 'exon':
            gtf_dict['type'].append(line_ls[2])
            gtf_dict['start'].append(int(line_ls[3]))
            gtf_dict['end'].append(int(line_ls[4]))
            gtf_dict['gene_name'].append(line.split('\t')[-1].strip().split(';')[3][12:-1])

    gtf_df = pd.DataFrame(gtf_dict)

    # get gene length
    gene_set = set(gtf_df.index)
    gene_length_dict = {'gene_name':[], 'length':[]}

    for gene in gene_set:
        intervals = list( np.array(gtf_df.loc[[gene]][['start','end']]) )
        merge_res = merge(intervals)
        length = sum([inval[1]-inval[0] for inval in merge_res])
        
        gene_length_dict['gene_name'].append(gene)
        gene_length_dict['length'].append(length)    
    gene_length_kb_df = pd.DataFrame(gene_length_dict).set_index('gene_name')/1e3

    # do normalization

    mapped_reads_Mb_Series = changhai_expr_mat.sum(axis=0) / 1e6

    useful_gene = set(gene_length_kb_df.index) & set(expr_mat.index)
    useful_expr_mat = expr_mat.loc[useful_gene]

    for gene in useful_expr_mat.index:
        for sample in useful_expr_mat.columns:
            useful_expr_mat.loc[gene, sample] = \
                useful_expr_mat.loc[gene, sample]/(gene_length_kb_df['length'][gene]*mapped_reads_Mb_Series[sample])

    log2_rpkm_changhai_expr_mat = np.log2(useful_expr_mat+0.1)

    return log2_rpkm_changhai_expr_mat


def remove_dup_genes(expr_mat, dup='first'):
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
    
    elif dup == 'first':  
        ## first
        unduplicated_expr_mat = expr_mat[~expr_mat.index.duplicated(keep='first')]

    elif dup == 'discard':  
        ## discard
        unduplicated_expr_mat = expr_mat[~dup_row]
        
    return unduplicated_expr_mat
    # return expr_mat[~expr_mat.index.duplicated(keep='first')]


def file_remove_dup_genes(exp_path, outpath, dup='first'):
    '''

    '''
    exp = pd.read_csv(exp_path, sep="\t", index_col=0, keep_default_na=False)
    exp_nodup = remove_dup_genes(expr_mat=exp, dup=dup)
    exp_nodup.to_csv(outpath, sep='\t')
    print('ok !!')

def remove_expr_less_genes():
    pass




def make_name_R_useful(name_array, symb_be_change='-', front_add = 'X'):
    new_name_array = [front_add+item.replace(symb_be_change,'_') for item in name_array]
    return new_name_array