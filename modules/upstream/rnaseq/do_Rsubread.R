library(Rsubread)
library(tidyverse)

# align ===============
fastq_dir = '~/raw_data/PRJNA726826/fastq/'
fq_ls = list.files(fastq_dir,pattern='fastq.gz$')
fq_base_name = sapply(strsplit(fq_ls, "\\."),'[[',1)
fq_ls = paste(fastq_dir,fq_ls,sep='/')
# ===
# fq_ls


align.stat.ls = c()
for(i in 2:length(fq_ls)){
  cat('## now is',i)
  align.stat <- align(index="hg38",readfile1=fq_ls[i],
                      output_file=paste(fq_base_name[i],'BAM',sep='.'),
                      phredOffset=33,nthreads=4)
  align.stat.ls = c(align.stat.ls, align.stat)
}


# featureCount ===============

# bam_file_ls = c()
bam_dir = '~/raw_data/PRJNA726826/do_subread/'
bam_ls = list.files(bam_dir,pattern='BAM$')
bam_base_name = sapply(strsplit(bam_ls, "\\."),'[[',1)
bam_ls = paste(bam_dir,bam_ls,sep='/')

gtf_file = '/home/lgh/reference/hl_reference/gtf/ensembl/Homo_sapiens.GRCh38.98.gtf.gz'

GTF.attrType = 'gene_id' # maybe "gene_name"
# ===
# bam_ls




fc_res = featureCounts(files=bam_ls, annot.ext=gtf_file,isGTFAnnotationFile=TRUE,
                       GTF.featureType="exon",GTF.attrType="gene_id",nthreads=4)


tmp_exp_mat_count = as.data.frame(fc_res$counts)
write_tsv(tmp_exp_mat_count,'PRJNA726826_count_exp_mat.tsv')


#############################
# DownStream ===============
#############################


# get sample info

SraRunTable = read_csv('~/raw_data/PRJNA726826/SraRunTable.txt')

anno_table = SraRunTable[,c('Run','Sample Name')]
colnames(anno_table) = c('Run','SampleName')

anno_table$'group' = sapply(strsplit(anno_table$'SampleName', "\\-"),function(x) paste(x[1:3],collapse='-'))

# get express matrix

exp_mat = read_tsv('../do_subread/PRJNA726826_count_exp_mat.tsv')

# change colnames

colnames(exp_mat) = sapply(strsplit(colnames(exp_mat), "\\."),'[[',1)
colnames(exp_mat) = c('gene_name', (anno_table%>%filter(Run%in%colnames(exp_mat)[2:10]))$SampleName)



res_ls = list()

deseq2_exp_mat = (exp_mat%>%column_to_rownames('gene_name'))[( anno_table%>%filter(group%in%c('L02-PAOA-0h','L02-PAOA-12h')) )$SampleName]
condition_table = as.factor(( anno_table%>%filter(group%in%c('L02-PAOA-0h','L02-PAOA-12h')) )$'group') 
# ===
# ii                                       
dds <- DESeqDataSetFromMatrix(deseq2_exp_mat, 
                              DataFrame(condition_table), 
                              design= ~ condition_table)
dds <- dds[rowSums(counts(dds)) > 1,]
dds2 <- DESeq(dds)
resultsNames(dds2)
# acquire the results using function results(), and assign to res
res <- results(dds2)

summary(res)
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)

res_ls$X_0_12 = resdata






