
dir_ls = list.dirs(path ='/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/AML_omic/wes/tumor/',
                    full.names = TRUE, recursive = FALSE)

fq_path_ls = c()

for (dir in dir_ls){
    fq_ls = list.files(dir,pattern='fq.gz$')
    sample_ls = unique( sapply(lapply(fq_ls,function(x){ substring(x,1,nchar(x)-8)}) , '[[', 1) )
    #print(fq_ls)
    #print(sample_ls)
    for(samp in sample_ls){
        fq_path_ls = c(fq_path_ls, 
                       paste(c( dir,'/',samp,'_1.fq.gz' ),collapse=''),
                       paste(c( dir,'/',samp,'_2.fq.gz' ),collapse=''))
        
        if(file.exists(paste(c( dir,'/',samp,'_1.fq.gz' ),collapse=''))==FALSE){
            print('FALSE')
        }
        if(file.exists(paste(c( dir,'/',samp,'_2.fq.gz' ),collapse=''))==FALSE){
            print('FALSE')
        }    
    
    }
    # break
}

fq_path_ls





































