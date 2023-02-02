# 模板 http://192.168.137.2:6424/tree/my_home/bigdisk/workspace/thh



dir = "~/bigdisk/thh/MJ20211104148-FX2021112300138-涂鸿浩-中国科学院上海药物研究所-真核转录组纯测序-12个样品-20211201/"
ddir = dir

file_ls = list.files(dir,'.raw.fastq.gz$')
file_path_ls = file.path(dir,file_ls)

sample_names = as.vector(sapply(file_ls,function(x) substr(x,1,nchar(x)-16)))


sample_names = unique(sample_names)
cat(paste(sample_names,': ','--',sep=''), sep='\n')




### ===============
### POST ANALYSIS
### ===============


get_sample_names = function(df){
    
    colnames(df)[2:length(colnames(df))] = 
        sapply(strsplit(colnames(df)[2:length(colnames(df))],split='/'),'[[',2)
    colnames(df)[1] = 'Hugo_Symbol'
    return( df )
    
}


result_dir = '/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/raw_data/wwb/MJ20220307202-FX2022032300143-吴文彪-中国科学院上海药物研究所-真核转录组纯测序-12个样品-20220323/analysis/result/fc.txt'            

count_df = read_tsv(result_dir,skip=1)
count_df = count_df[,-c(2:6)] # 去掉featureCount前面的列

count_df = get_sample_names(count_df)
# count_df <- aggregate(.~ Hugo_Symbol, count_df, mean) # 解决基因名不唯一


deseq2_exp_mat = count_df[,c(1,2:4,5:7)]%>%column_to_rownames('Hugo_Symbol')
group_df = data.frame(
    sample=colnames(deseq2_exp_mat)[-1],
    group=c(rep('Vec',3),rep('GLS',3))
)
deg_gls_vec_Pair = 
    do_deseq2_new(deseq2_exp_mat,group_df,g1='Vec',g2='GLS')



deseq2_exp_mat = count_df[,c(1,8:10,11:13)]%>%column_to_rownames('Hugo_Symbol')
group_df = data.frame(
    sample=colnames(deseq2_exp_mat)[-1],
    group=c(rep('Vec',3),rep('GLS',3))
)
deg_gls_vec_EtOH = 
    do_deseq2_new(deseq2_exp_mat,group_df,g1='Vec',g2='GLS')



deg_gls_vec_EtOH$resdata

