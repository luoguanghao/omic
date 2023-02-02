import numpy as np
import pandas as pd
from pathlib import Path
import os, sys



def get_fq_file(fastq_path):
    # input fq file
    
    all_file = pd.DataFrame(os.walk(fastq_path))[2][0]

    # get file
    single_end_ls = []
    pair_end_ls = []
    for file in all_file:
        if 'fastq' in file:
            if '_' in file:
                pair_end_ls.append(file)
            else:
                single_end_ls.append(file)

    pair_end_ls = sorted(pair_end_ls)
    single_end_ls = sorted(single_end_ls)
    
    return {'pair_end_ls':pair_end_ls, 'single_end_ls':single_end_ls}


def read_excel():
    '''
    如何pandas读取excel
    https://www.jianshu.com/p/c3c2ac84fb02
    '''
    pass


def read_txt(file_path, sep='\t', LF='\n', columns_name_line=0, file_type=None):
    '''
    columns_name_line: which line to be set as columns name, start from 0; the next line of columns name would be set as start of data
    '''

    file = open(file_path).read().strip().split('%s'%LF)
    file = [ line.split('%s'%sep) for line in file ]

    if file_type == 'GDSC':
        for i in range(len(file)):
            if len(file[i])==14:
                file[i] = file[i][0:1]+file[i][2:]
            elif len(file[i])==15:
                file[i] = file[i][0:1]+file[i][3:]

    columns_name = file[columns_name_line]
    
    df = pd.DataFrame(file[columns_name_line+1:], columns=columns_name)

    return df



def do_wget():
    pass

def do_curl():
    pass












