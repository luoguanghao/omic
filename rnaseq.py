from pathlib import Path
import os, subprocess
import numpy as np
import pandas as pd



def do_one_process(SRR_no, fq_1_path, fq_2_path, reference_genome_path, gtf_path, work_dir):
    
    os.chdir(work_dir)
    
    ##########
    # HISAT2 #
    ##########
    '''
    hisat2 -t -p 10 -x /y/Bondi/work/try/try1/hg38/genome \
    -1 ./SRR5381445_1.fastq -2 ./SRR5381445_2.fastq -S ./SRR5381445.sam
    '''

    os.mkdir('map')

    now_dir = '%s%s/'%(work_dir,'map')
    sam_path = now_dir+'%s.sam'%SRR_no  #########
    log_path = now_dir+'hisat2.log'

    cmd = 'hisat2 -t -p 10 -x %s -1 %s -2 %s -S %s > %s 2>&1'%(reference_genome_path, fq_1_path, fq_2_path, sam_path, log_path)

    p = subprocess.Popen(cmd, shell=True)

    # print(p.stdout.read())
    # output = p.communicate()
    # print(output)

    status = p.wait()

    print('done!!\n'+cmd)
    ##########
    ## sort ##
    ##########
    os.mkdir('sort')

    now_dir = '%s%s/'%(work_dir,'sort')
    bam_path = now_dir+'%s.bam'%SRR_no
    sort_bam_path = now_dir+'%s.sorted.bam'%SRR_no ##########
    log_path = now_dir+'error.log'

    cmd = 'samtools view -S %s -b > %s'%(sam_path, bam_path)

    p = subprocess.Popen(cmd, shell=True)
    status = p.wait()

    print('done!!\n'+cmd)


    cmd = 'samtools sort %s -o %s 2>> %s'%(bam_path, sort_bam_path, log_path) ####

    p = subprocess.Popen(cmd, shell=True)
    status = p.wait()

    print('done!!\n'+cmd)

    ###########
    ## index ##
    ###########
    cmd = 'samtools index %s 2>> %s'%(sort_bam_path, log_path) ####

    p = subprocess.Popen(cmd, shell=True)
    status = p.wait()

    print('done!!\n'+cmd)

    ###############
    # HTSeq-count #
    ###############
    os.mkdir('count')

    now_dir = '%s%s/'%(work_dir,'count')
    exp_count_path = now_dir+'%s_exp_htcount.tsv'%SRR_no
    log_path = now_dir+'error.log'

    # cmd = 'htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id -m union %s %s > %s 2> %s'%(sort_bam_path, gtf_path, exp_count_path, log_path)
    cmd = 'htseq-count -f bam -r name -s no -i gene_name %s %s > %s 2> %s'%(sort_bam_path, gtf_path, exp_count_path, log_path)

    p = subprocess.Popen(cmd, shell=True)

    status = p.wait()

    print('done!!\n'+cmd)

    ###### ^-^ ######
    print('All Done!!!^-^')


    
def multi_run_process():
    pass

    
def combine_into_exp_mat():
    pass




if __name__ == '__main__':
    
    ####################################
    #  20个20个的跑
    ####################################
    
    work_dir = '/y/Bondi/work/try/try3/'
    # os.chdir(work_dir)

    hisat2_path = '/y/home/lgh/bin/hisat2'
    samtools_path = '/y/home/lgh/bin/samtools'
    htseq_count_path = '/y/home/lgh/.local/bin/htseq-count'

    reference_genome_path = '/y/home/lgh/reference/index/hisat/hg38/genome'
    gtf_path = '/y/home/lgh/reference/gtf_file/Homo_sapiens.GRCh38.89.chr.gtf'

    SRR_no = 'SRR5381447'
    fq_1_path = '/y/Bondi/data/Qiu_Hui_2019_LIMORE/GSE97098/fastq_file/SRR5381447_1.fastq'
    fq_2_path = '/y/Bondi/data/Qiu_Hui_2019_LIMORE/GSE97098/fastq_file/SRR5381447_2.fastq'

    
    do_one_process(SRR_no, fq_1_path, fq_2_path, reference_genome_path, gtf_path, work_dir)






    