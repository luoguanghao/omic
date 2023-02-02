library(filesstrings)


fastq_dir = '/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/AML_omic/wes/tumor/'
fq_ls = list.files(fastq_dir,pattern='fq.gz$')

sample_ls = unique( sapply(lapply(fq_ls,function(x){ substring(x,1,nchar(x)-8)}) , '[[', 1) )

wd = paste(getwd(),'tumor/',sep='/')

# mv
for (i in 1:length(sample_ls)){
    new_dir = paste(wd,sample_ls[i],sep='/')

    dir.create(new_dir)
    file_name = paste(sample_ls[i],'_1.fq.gz',sep='')
    file_path = paste(wd,file_name,sep='')
    file.move(file_path, new_dir) # mv from file_path to new_dir
    cat(paste(c('@@ move',file_path,'to',new_dir),collapse=' '))
    file_name = paste(sample_ls[i],'_2.fq.gz',sep='')
    file_path = paste(wd,file_name,sep='')
    file.move(file_path, new_dir)
    cat(paste(c('@@ move',file_path,'to',new_dir),collapse=' '))
    }















