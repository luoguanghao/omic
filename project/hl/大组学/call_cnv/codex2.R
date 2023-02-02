
library(CODEX2)
library(WES.1KG.WUGSC) # Load Toy data from the 1000 Genomes Project.


dirPath <- '/home/lgh/my_project/hl_call_cnv/data/pdcmtation-2021-3-8批次_2/bam'

bamFile <- list.files(dirPath, pattern ='*.bam$')

bamdir <- file.path(dirPath, bamFile)
# get bam file
## ctrl
dirPath_J012_ctrl <- '/media/bio/汪翰林硬盘二代目/2021-AML-omics/wes-ctrl/X101SC20061552-Z01-J012/02.Bam'
bamFile_J012_ctrl <- list.files(dirPath_J012_ctrl, pattern ='*.bam$')
bamdir_J012_ctrl <- file.path(dirPath_J012_ctrl, bamFile_J012_ctrl)


dirPath_J011_ctrl <- '/media/bio/汪翰林硬盘二代目/2021-AML-omics/wes-ctrl/ctrl_MUTATION-0513/X101SC20061552-Z01-J011/02.Bam'
bamFile_J011_ctrl <- list.files(dirPath_J011_ctrl, pattern ='*.bam$')
bamdir_J011_ctrl <- file.path(dirPath_J011_ctrl, bamFile_J011_ctrl)


dirPath_J010_ctrl <- '/media/bio/汪翰林硬盘二代目/2021-AML-omics/wes-ctrl/ctrl_mutation-amlhealthy-4-27/X101SC20061552-Z01-J010/02.Bam'
bamFile_J010_ctrl <- list.files(dirPath_J010_ctrl, pattern ='*.bam$')
bamdir_J010_ctrl <- file.path(dirPath_J010_ctrl, bamFile_J010_ctrl)

#head(bamFile_J010_ctrl)

#head(bamFile_J011_ctrl)

#head(bamFile_J012_ctrl)

sample_names = sapply(strsplit(multiannoFile,'\\.'),'[[',1)


## tumor
dirPath_re_tm <- '/media/bio/汪翰林硬盘二代目/2021-AML-omics/wes-tumor/tumor_aml-PDC-mutation-2021-4-1 補充/re/02.Bam'
bamFile_re_tm <- list.files(dirPath_re_tm, pattern ='*.bam$')
bamdir_re_tm <- file.path(dirPath_re_tm, bamFile_re_tm)


dirPath_J006_tm <- '/media/bio/汪翰林硬盘二代目/2021-AML-omics/wes-tumor/tumor_aml-pdc-mutation-2021-4-4/liuchang0401026/Results-X101SC20061552-Z01-J006-B52-27-20210329/02.Bam'
bamFile_J006_tm <- list.files(dirPath_J006_tm, pattern ='*.bam$')
bamdir_J006_tm <- file.path(dirPath_J006_tm, bamFile_J006_tm)


dirPath_J008_tm <- '/media/bio/汪翰林硬盘二代目/2021-AML-omics/wes-tumor/tumor_mutation-4-19/008/results/02.Bam'
bamFile_J008_tm <- list.files(dirPath_J008_tm, pattern ='*.bam$')
bamdir_J008_tm <- file.path(dirPath_J008_tm, bamFile_J008_tm)


dirPath_J009_tm <- '/media/bio/汪翰林硬盘二代目/2021-AML-omics/wes-tumor/tumor_mutation-4-25/X101SC20061552-Z01-J009/02.Bam'
bamFile_J009_tm <- list.files(dirPath_J009_tm, pattern ='*.bam$')
bamdir_J009_tm <- file.path(dirPath_J009_tm, bamFile_J009_tm)


## bed file
bedFile <- file.path('~/my_upstream_project/hl_call_cnv/data/novo_wes.bed')


## combine data

bamFile = c(bamFile_re_tm,bamFile_J006_tm,bamFile_J008_tm,bamFile_J009_tm,
            bamFile_J012_ctrl,bamFile_J011_ctrl,bamFile_J010_ctrl)
bamdir = c(bamdir_re_tm,bamdir_J006_tm,bamdir_J008_tm,bamdir_J009_tm,
            bamdir_J012_ctrl,bamdir_J011_ctrl,bamdir_J010_ctrl)

sampname <- sapply(strsplit(bamFile, "\\."),'[[',1)


# getbambed
bambedObj <- getbambed(bamdir = bamdir, bedFile = bedFile, 
                       sampname = sampname, projectname = "CODEX2_demo")  #<<<<<<<<<<<<<<<<<<<<
bamdir <- bambedObj$bamdir; sampname <- bambedObj$sampname
ref <- bambedObj$ref; projectname <- bambedObj$projectname

head(ref)
str(ref)


## Getting GC content and mappability ##
genome = BSgenome.Hsapiens.UCSC.hg19 # hg19
# library(BSgenome.Hsapiens.UCSC.hg38); genome = BSgenome.Hsapiens.UCSC.hg38 # hg38
gc <- getgc(ref, genome = genome)
mapp <- getmapp(ref, genome = genome)
values(ref) <- cbind(values(ref), DataFrame(gc, mapp))  


## Getting raw read depth ##
coverageObj <- getcoverage(bambedObj, mapqthres = 20)
Y <- coverageObj$Y
write.csv(Y, file = paste(projectname, '_coverage.csv', sep=''), quote = FALSE)


head(Y)


## Quality control ##
qcObj <- qc(Y, sampname, ref, cov_thresh = c(20, 4000),
            length_thresh = c(20, 2000), mapp_thresh = 0.9,
            gc_thresh = c(20, 80))
Y_qc <- qcObj$Y_qc; sampname_qc <- qcObj$sampname_qc
ref_qc <- qcObj$ref_qc; qcmat <- qcObj$qcmat; gc_qc <- ref_qc$gc
write.table(qcmat, file = paste(projectname, '_qcmat', '.txt', sep=''),
            sep = '\t', quote = FALSE, row.names = FALSE)


head(ref)


####################
## Running CODEX2 ##
####################

colnames(Y_qc) <- paste('sample_', 1:ncol(Y_qc), sep='') # <<<


Y.nonzero <- Y_qc[apply(Y_qc, 1, function(x){!any(x==0)}),]
pseudo.sample <- apply(Y.nonzero,1,function(x){exp(1/length(x)*sum(log(x)))})
N <- apply(apply(Y.nonzero, 2, function(x){x/pseudo.sample}), 2, median)



## Running CODEX2 with negative control samples ##
chr <- 20 # < < < < < <
chr.index <- which(seqnames(ref_qc)==chr)
# Below are pre-computed demo dataset, stored as part of the CODEX2 R-package.
normObj <- normalize_codex2_ns(Y_qc = Y_qc[chr.index,],
                               gc_qc = gc_qc[chr.index], 
                               K = 1:5, norm_index = norm_index_demo,  # 正常样本就是在这里引入
                               N = N)

norm_index_demo



Yhat.ns <- normObj$Yhat; fGC.hat.ns <- normObj$fGC.hat;
beta.hat.ns <- normObj$beta.hat; g.hat.ns <- normObj$g.hat; h.hat.ns <- normObj$h.hat
AIC.ns <- normObj$AIC; BIC.ns <- normObj$BIC; RSS.ns <- normObj$RSS


choiceofK(AIC.ns, BIC.ns, RSS.ns, K = 1:5 , filename = "codex2_ns_choiceofK.pdf")






# Running segmentation by CODEX2
finalcall.CBS <- segmentCBS(Y_qc[chr.index,],  # recommended
                            Yhat.ns, optK = which.max(BIC.ns),
                            K = 1:5,
                            sampname_qc = colnames(Y_qc),
                            ref_qc = ranges(ref_qc)[chr.index],
                            chr = chr, lmax = 400, mode = "integer")

#finalcall.HMM <- segmentHMM(Y_qc[chr.index,],  # not recommended
#                            Yhat.ns, optK = which.max(BIC.ns),
#                            K = 1:5,
#                            sampname_qc = colnames(Y_qc),
#                            ref_qc = ranges(ref_qc)[chr.index],
#                            chr = chr, mode = "integer")


# Post-segmentation pruning and filtering are recommended based on CNV length (filter1), length per exon (filter2), likelihood ratio (filter3), and number of exons (filter4).

filter1 <- finalcall.CBS$length_kb<=200
filter2 <- finalcall.CBS$length_kb/(finalcall.CBS$ed_exon-finalcall.CBS$st_exon+1)<50
finalcall.CBS.filter <- finalcall.CBS[filter1 & filter2, ]

filter3 <- finalcall.CBS.filter$lratio>40
filter4 <- (finalcall.CBS.filter$ed_exon-finalcall.CBS.filter$st_exon)>1
finalcall.CBS.filter_final=finalcall.CBS.filter_final[filter3|filter4,]




View(finalcall.CBS.filter_final)




























