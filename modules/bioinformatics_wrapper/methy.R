

##
# group_df
# methy_path = "~/my_project/data/TCGA/BRCA/TCGA-BRCA.methylation450.tsv.gz"
##

load_TCGA_methy_and_preprocess <- function(group_df, methy_path){
    library(data.table)
    library(impute)
    library(ChAMP)
    library(stringr)
    library(tibble)

    methy_df = data.table::fread(methy_path, data.table = F)
    methy_df = column_to_rownames(methy_df,"Composite Element REF")
    colnames(methy_df)= str_sub(colnames(methy_df),1,15)

    ## preprocess , 补NA值
    beta=as.matrix(methy_df)
    # beta信号值矩阵里面不能有NA值
    beta=impute.knn(beta) 
    sum(is.na(beta))

    beta=beta$data
    beta=beta+0.00001

    ## load into champ
    pd = data.frame(Sample_Name=colnames(beta))
    pd = merge(pd,group_df,all.x=TRUE,by.x='Sample_Name',by.y='sample_id_B')
    pd = pd[!duplicated(pd$Sample_Name),]
    pd$group[is.na(pd$group)]='N'
    pd$group = as.factor(pd$group)

    myLoad=champ.filter(beta = beta[,pd$Sample_Name], pd=pd) #这一步已经自动完成了过滤

    # normalization
    myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=8)
    # 归一化过程产生了缺失值,需要将有NA的样本和它们的配对样本一起删掉
    num.na <- apply(myNorm,2,function(x)(sum(is.na(x))))
    table(num.na)

    ## 把有na的样本去掉
    names(num.na) = colnames(myNorm)
    dt = names(num.na[num.na>0])
    #dn = str_replace(dt,"-01","-11")
    keep = setdiff(colnames(myNorm),dt)
    myNorm = myNorm[,keep]
    pd = myLoad$pd
    pd = pd[pd$Sample_Name %in% keep,]
    identical(pd$Sample_Name,colnames(myNorm))

    return(list(pd=pd, myNorm=myNorm, beta=beta))
} 









