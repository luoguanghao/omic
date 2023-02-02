# 文件处理
#
# https://blog.csdn.net/weixin_40628687/article/details/79249861
# 
# 
##  
## 
# 输入:
# fastq_dir fq文件目录
# wd 工作目录
# 这里就是把fastq_dir里面的fq文件根据其样本归属，分别在wd建立文件夹，然后把fq文件移进去
fastq_dir = '/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/AML_omic/wes/tumor/'
wd = paste(getwd(),'tumor/',sep='/')


# ===============
library(filesstrings)


fq_ls = list.files(fastq_dir,pattern='fq.gz$')

sample_ls = unique( sapply(lapply(fq_ls,function(x){ substring(x,1,nchar(x)-8)}) , '[[', 1) )


for (i in 1:length(sample_ls)){
    new_dir = paste(wd,sample_ls[i],sep='/')
    dir.create(new_dir)
    file_name = paste(sample_ls[i],'_1.fq.gz',sep='')
    file_path = paste(wd,file_name,sep='')
    file.move(file_path, new_dir)
    cat(paste(c('@@ move',file_path,'to',new_dir),collapse=' '))
    file_name = paste(sample_ls[i],'_2.fq.gz',sep='')
    file_path = paste(wd,file_name,sep='')
    file.move(file_path, new_dir)
    cat(paste(c('@@ move',file_path,'to',new_dir),collapse=' '))
}



































'/my_upstream_project/hl_call_cnv/data/fq_file/tumor/DCY20181115'





