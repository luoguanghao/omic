# 找到所有tumor bam,输出yaml需要的格式: sample_name, dir

## tumor
dirPath_re_tm <- '/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/AML_omic/wes/tumor/bam'
bamFile_re_tm <- list.files(dirPath_re_tm, pattern ='*.bam$')
bamdir_re_tm <- file.path(dirPath_re_tm, bamFile_re_tm)

sample_ls = sapply(strsplit(bamFile_re_tm,'\\.'),'[[',1)

bamdir_re_tm_dockerpath = paste('/gatk/top/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/AML_omic/wes/tumor/bam/',bamFile_re_tm,sep='')

cat(paste( paste(sample_ls, bamdir_re_tm_dockerpath, sep=': '), collapse='\n' ))













