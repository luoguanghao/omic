

## PE

fq_dir1 = '/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/YQ//MJ20210902153-FX2021091600035-余强-中国科学院上海药物研究所-医学有参转录组测序-9个样品-20210916/'  
fq_ls1 = list.files(fq_dir1,'gz$') # 获取某个后缀的文件


index = '/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/raw_data/rnaseq_genome/grcm38/genome'
ref = '-'
gtf = '/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/raw_data/rnaseq_genome/Mus_musculus.GRCm38.94.chr_patch_hapl_scaff.gtf'

#sample = ''
#fq1 = ''
#fq2 = ''
samples = unique(as.vector(sapply(fq_ls1,function(x) substr(x,1,nchar(x)-16))))

for(s in samples){
    fq1 = paste(fq_dir1,s,'.R1.raw.fastq.gz',sep='')
    fq2 = paste(fq_dir1,s,'.R2.raw.fastq.gz',sep='')
    log = sprintf('%s.log',s)
    cat(paste('nohup bash /home/lgh/my_project/omic/modules/upstream/rnaseq/hisat2_samtools_pe.sh', 
        '-i',index , '-r',ref , '-g',gtf , '-n',s , 
        '-a',fq1 , '-b',fq2,'>', log,'&'))
    cat('\n\n')
}


# =========


# 批量FeatureCount
for(s in samples){
    bam = paste(s,'/',s,'_sorted.bam',sep='')
    output = paste(s,'/',s,'_fc.txt',sep='')
    log = sprintf('%s_fc.log',s)
    cat(paste('nohup featureCounts', '-T 4 -p -O -g gene_name', 
            '-a',gtf , '-o',output, bam, '>', log,'&'))
    cat('\n\n')
}

# ======== featureCount all ======
fq_ls1 = ''
fq_ls2 = ''

samples1 = unique(as.vector(sapply(fq_ls1,function(x) substr(x,1,nchar(x)-16))))
samples2 = unique(as.vector(sapply(fq_ls2,function(x) substr(x,1,nchar(x)-16))))                                  

bam_ls = c()
for(s in samples1){
bam_ls = c(bam_ls,paste('/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/YQ/upstream_20210916/',s,'/',s,'_sorted.bam',sep=''))
}
for(s in samples2){
bam_ls = c(bam_ls,paste('/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/YQ/upstream_20210922/',s,'/',s,'_sorted.bam',sep=''))
}
                                   
output = paste('fc.txt')
log = sprintf('fc.log',s)
cat(paste('featureCounts', '-T 8 -p -O -g gene_name', '-a',gtf , '-o',output, paste(bam_ls,collapse=' ')))   


# ==== 修饰featureCount结果 ====
df = read_tsv('/media/bio/89c80d9d-94c7-47e2-89b8-515d83753ae5/server/YQ/featurecount/fc.txt',skip=1)  
colnames(df)[7:length(colnames(df))] = sapply(strsplit(colnames(df)[7:length(colnames(df))],'/'),'[[',8)

write_tsv(df,'../featurecount/fc_yq.tsv')








