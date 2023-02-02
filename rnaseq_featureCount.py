'''
problem:
bin path,没有给
函数的返回值 do函数

没有给函数写说明
'''
import sys
sys.path.append("/y/Bondi/src/")
from utilities import file_tools, NCBI_tools

from pathlib import Path
import os, subprocess
import numpy as np
import pandas as pd
import multiprocessing
import time

def new_hisat2_core(reference_genome_path, fq_path, sam_path, log_path, new_folder=True, new_dir=None, t=4):

    print('[Action!]')
    
    if new_folder == True:
        os.mkdir(new_dir)
    
    if str(type(fq_path)) == "<class 'list'>":
        cmd = 'hisat2 -t -p %s -x %s -1 %s -2 %s -S %s > %s 2>&1'%(t, reference_genome_path, fq_path[0], fq_path[1], sam_path, log_path)
    else:
        cmd = 'hisat2 -t -p %s -x %s -U %s -S %s > %s 2>&1'%(t, reference_genome_path, fq_path, sam_path, log_path)

    p = subprocess.Popen(cmd, shell=True)

    status = p.wait()
    # time.sleep(3)
    
    return {'status':1, 'cmd':cmd}


def hisat2_core(single_pair, reference_genome_path, log_path, sam_path, fq_1_path=None, fq_2_path=None, fq_path=None, t=4):

    if single_pair == 'pair':
        cmd = 'hisat2 -t -p %s -x %s -1 %s -2 %s -S %s > %s 2>&1'%(t, reference_genome_path, fq_1_path, fq_2_path, sam_path, log_path)
    else:
        cmd = 'hisat2 -t -p %s -x %s -U %s -S %s > %s 2>&1'%(t, reference_genome_path, fq_path, sam_path, log_path)



    p = subprocess.Popen(cmd, shell=True)

    status = p.wait()

    return {'status':status, 'cmd':cmd}
    

def do_hisat2(fastq_name, single_pair, work_dir, hisat2_path, reference_genome_path, t):
    
    if single_pair == 'pair':
        SRR_no = fastq_name.split('_')[0]
        
    else:
        SRR_no = fastq_name.split('.')[0]
    # SRR_no = 'SRR5381447'
    
    print('[Action!] start: %s type is %s...'%(SRR_no, single_pair))
    
    now_dir = '%s/%s/'%(work_dir, SRR_no)
    
    # print('### ',now_dir)
    
    # ======== #
    os.mkdir(now_dir)
    # print('[New!] Create folder %s'%now_dir)
    os.chdir(now_dir)
   
    

    ##########
    # HISAT2 #
    ##########
    
    # prepare for hisat2
    if single_pair == 'pair':
        fq_1_path = '%s%s_1.fastq'%(fastq_path,SRR_no)
        fq_2_path = '%s%s_2.fastq'%(fastq_path,SRR_no)
    else:
        fq_path = '%s%s.fastq'%(fastq_path,SRR_no)

    sam_path = now_dir+'%s.sam'%SRR_no  #########
    hisat2_log_path = now_dir+'%s_hisat2.log'%SRR_no
    
    
    
    # action!
    if single_pair == 'pair':
        do_hisat2_res = hisat2_core(single_pair=single_pair, log_path=hisat2_log_path, 
                         reference_genome_path=reference_genome_path, sam_path=sam_path, 
                         fq_1_path=fq_1_path, fq_2_path=fq_2_path, t=t)
    else:
        do_hisat2_res = hisat2_core(single_pair=single_pair, log_path=hisat2_log_path, 
                         reference_genome_path=reference_genome_path, sam_path=sam_path, 
                         fq_path=fq_path, t=t)        
    
    # print('[Return] %s'%do_res['status'])
    print('[done]\n'+do_hisat2_res['cmd']+'\n[done]')


def new_multi_run_hisat2(single_end_fq_ls, pair_end_fq_ls, fastq_path, work_dir, 
                             reference_genome_path, hisat2_path=None,t=4 , p=3):
    SRR_no_ls = []
    new_dir_ls = []
    fq_path_ls = []
    sam_path_ls = []
    hisat2_log_path_ls = []
    
    for i_fq_file in range(0, len(pair_end_fq_ls), 2):
        SRR_no = pair_end_fq_ls[i_fq_file].split('_')[0]
        SRR_no_ls.append(SRR_no)
        new_dir = '%s/%s/'%(work_dir,SRR_no)
        new_dir_ls.append(new_dir)
        
        sam_path_ls.append(new_dir+'%s.sam'%SRR_no)
        fq_path_ls.append( ['%s%s_1.fastq'%(fastq_path,SRR_no), '%s%s_2.fastq'%(fastq_path,SRR_no)] )
        hisat2_log_path_ls.append(new_dir+'%s_hisat2.log'%SRR_no)
        
    for i_fq_file in range(0,len(single_end_fq_ls), 1):
        SRR_no = single_end_fq_ls[i_fq_file].split('.')[0]
        SRR_no_ls.append(SRR_no)
        new_dir = '%s/%s/'%(work_dir,SRR_no)
        new_dir_ls.append(new_dir)
        
        sam_path_ls.append(new_dir+'%s.sam'%SRR_no)
        fq_path_ls.append( '%s%s_1.fastq'%(fastq_path,SRR_no) )
        hisat2_log_path_ls.append(new_dir+'%s_hisat2.log'%SRR_no)
    
    # print(new_dir_ls)
    
    # multi run
    pool = multiprocessing.Pool(processes = p)
    res_ls = []        
    for i in range(len(fq_path_ls)):
        # do_fc_res = featurecount_core(gtf_path, count_path_ls[i], sam_path_ls[i], fc_log_path_ls[i])
        kwargs = {
                   'reference_genome_path':reference_genome_path, 'fq_path':fq_path_ls[i], 
                'sam_path':sam_path_ls[i], 'log_path':hisat2_log_path_ls[i], 'new_dir':new_dir_ls[i], 't':t
                }
        res = pool.apply_async(new_hisat2_core, args=(), kwds=kwargs)
        res_ls.append(res)
    
    pool.close()
    pool.join() # 调用join之前，先调用close函数，否则会出错。
    
    return res_ls        
                 

def multi_run_hi(fastq_ls, single_pair, work_dir, hisat2_path, reference_genome_path, t, p):
    
    print('$$$ Multi-Run : processes=%s $$$'%p)
    pool = multiprocessing.Pool(processes = p)
    
    if single_pair is 'pair':

        for i_fn in range(0, len(fastq_ls),2):
            
            #print('... %s')%i_fn
            
            fastq_name = fastq_ls[i_fn]
            kwargs = {'fastq_name':None, 'single_pair':'pair', 'work_dir':work_dir, 
                'hisat2_path':hisat2_path, 'reference_genome_path':reference_genome_path, 't':t}       
            kwargs['fastq_name'] = fastq_name
    
            
    
            #pool.apply_async(func, (fastq_name,))
            pool.apply_async(do_hisat2, args=(), kwds=kwargs)

        pool.close()
        pool.join() # 调用join之前，先调用close函数，否则会出错。
    else:
        
        for i_fn in range(0, len(fastq_ls)):

            fastq_name = fastq_ls[i_fn]
            kwargs = {'fastq_name':None, 'single_pair':'single', 'work_dir':work_dir, 
                'hisat2_path':hisat2_path, 'reference_genome_path':reference_genome_path, 't':t}       
            kwargs['fastq_name'] = fastq_name

            #pool.apply_async(func, (fastq_name,))
            pool.apply_async(do_hisat2, args=(), kwds=kwargs)

        pool.close()
        pool.join() # 调用join之前，先调用close函数，否则会出错。        


def featurecount_core(gtf_path, count_path, sam_path, log_path):

    print('[Action!]')
    
    cmd = 'featureCounts -p -t exon -g gene_id -T 10 -a %s -o %s %s > %s 2>&1'%(gtf_path, count_path, sam_path, log_path)

    p = subprocess.Popen(cmd, shell=True)


    status = p.wait()
    
    # time.sleep(3)
    print('[Done]\n%s\n[Done]\n'%cmd)
    return {'status':status, 'cmd':cmd}


def multi_run_featurecount(single_end_fq_ls, pair_end_fq_ls, work_dir, gtf_path, featureCount_path=None, p=3):

    sam_path_ls = []
    count_path_ls = []
    fc_log_path_ls = []
    for i_fq_file in range(0, len(pair_end_fq_ls), 2):
        SRR_no = pair_end_fq_ls[i_fq_file].split('_')[0]
        now_dir = '%s/%s/'%(work_dir,SRR_no)
        sam_path_ls.append(now_dir+'%s.sam'%SRR_no)
        count_path_ls.append(now_dir+'%s_featureCount.tsv'%SRR_no)
        fc_log_path_ls.append(now_dir+'%s_featureCount.log'%SRR_no)
    for i_fq_file in range(0,len(single_end_fq_ls), 1):
        SRR_no = single_end_fq_ls[i_fq_file].split('.')[0]
        now_dir = '%s/%s/'%(work_dir,SRR_no)
        sam_path_ls.append(now_dir+'%s.sam'%SRR_no)
        count_path_ls.append(now_dir+'%s_featureCount.tsv'%SRR_no)
        fc_log_path_ls.append(now_dir+'%s_featureCount.log'%SRR_no)    
    
    # multi run
    pool = multiprocessing.Pool(processes = p)
    res_ls = []
    for i in range(len(sam_path_ls)):
        # do_fc_res = featurecount_core(gtf_path, count_path_ls[i], sam_path_ls[i], fc_log_path_ls[i])
        kwargs = {
                'gtf_path':gtf_path, 'count_path':count_path_ls[i], 
                'sam_path':sam_path_ls[i], 'log_path':fc_log_path_ls[i]
                }
        res = pool.apply_async(featurecount_core, args=(), kwds=kwargs)
        res_ls.append(res)
    
    pool.close()
    pool.join() # 调用join之前，先调用close函数，否则会出错。
    
    return res_ls


def combine_expr_mat(work_dir, GSE_no, count_path_ls, SRR_no_ls, meta_data=None):
    count_df_ls = []

    for i_c_p in range(0,len(count_path_ls)):
       
        #sample_name = 
        counts_df = pd.read_csv(count_path_ls[i_c_p], sep='\t', header=1, index_col=0)
        
        counts_df_columns = list(counts_df.columns)
        if meta_data != None:
            counts_df_columns[-1] = meta_data.loc[SRR_no_ls[i_c_p],'sample_name']
        counts_df.columns = counts_df_columns
        
        count_df_ls.append(counts_df)  

    all_counts_df = pd.concat(
        [count_df_ls[0]]+[count_df_ls[i_cd].iloc[:,-1] for i_cd in range(1,len(count_df_ls))], axis=1, join='outer'
        )
    
    return all_counts_df


def get_path(GSE_no, work_dir, SRR_no, mk=False):
    pass


if __name__ == "__main__":

    # hasit2 important file

    reference_genome_path = '/y/home/lgh/reference/index/hisat/hg38/genome'
    gtf_path = '/y/home/lgh/reference/gtf_file/Homo_sapiens.GRCh38.89.chr.gtf'

    fastq_path = '/y/Bondi/data/Qiu_Hui_2019_LIMORE/GSE78236/fastq_file/'
    
    # get file
    get_fq_file_res = file_tools.get_fq_file(fastq_path)
    
    pair_end_ls = get_fq_file_res['pair_end_ls']
    single_end_ls = get_fq_file_res['single_end_ls']

    work_dir = '/y/Bondi/work/RNAseq/GSE97098_tmp/'

    # Action!
    res = new_multi_run_hisat2(single_end_ls, pair_end_ls, fastq_path, work_dir, reference_genome_path, t=8, p=5)

    multi_run_featurecount(single_end_ls, pair_end_ls, work_dir, gtf_path)



    