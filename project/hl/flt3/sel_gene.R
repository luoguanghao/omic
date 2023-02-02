# http://localhost:8888/notebooks/d/my_project/hl/vizome%E5%8F%91%E7%8E%B0/select_cebpa_score.ipynb
# 贪婪算法选基因

# source('/mnt/d/omic/modules/plot/survival.R')
source('/mnt/d/omic/project/hl/hl_分测试集验证集_func版.R')
library("openxlsx")
library(DESeq2)
library(data.table)

cal_corr <- function(exp_mat,drug,drug_sen_df,raw_cebpa_set){
    # drug
    # i_drug = 1
    x_mat = t(exp_mat[raw_cebpa_set,intersect(colnames(exp_mat), (drug_sen_df%>%filter(inhibitor==drug))$lab_id)])
    y_vec = (drug_sen_df%>%filter(inhibitor==drug)%>%column_to_rownames('lab_id'))[
                intersect(colnames(exp_mat), (drug_sen_df%>%filter(inhibitor==drug))$lab_id),]$auc

    corr_df = list(gene=c(),r=c(),p=c())
    for(i in 1:dim(x_mat)[2]){
        tmp = cor.test(x_mat[,i],y_vec)
        corr_df$gene = c(corr_df$gene,colnames(x_mat)[i])
        corr_df$r = c(corr_df$r,tmp$estimate)
        corr_df$p = c(corr_df$p,tmp$p.value)
    }
    corr_df = data.frame(corr_df)
    rownames(corr_df) = NULL
    Mido_corr_df = corr_df[order(corr_df$p),]

    Mido_corr_gene = Mido_corr_df$gene[1:50]
    return(list(Mido_corr_df=Mido_corr_df,Mido_corr_gene=Mido_corr_gene))
}

### prepare data
## mut data
maf_df = as.data.frame(fread('../../data/VIZOME/LAML_nature_PDC.maf',sep='\t'))
ITD_sample_ls = as.vector((maf_df %>% dplyr::filter(Hugo_Symbol=='FLT3'&ITDorNOT=='ITD'))['Tumor_Sample_Barcode'][['Tumor_Sample_Barcode']])
ITD_sample_ls = str_replace_all(paste('X',ITD_sample_ls,sep=''),'-','.')

## exp data
#exp_mat = as.data.frame(fread('../../data/VIZOME/nature_aml_log2_fpkm.txt',sep='\t'))
# 处理重复基因，滤出需要的基因
#exp_mat = aggregate(.~ Symbol, exp_mat, mean)
#exp_mat = exp_mat %>% column_to_rownames('Symbol')
#colnames(exp_mat) = str_replace_all(paste('X',colnames(exp_mat),sep=''),'-','.')
exp_mat = read.table('../../data/VIZOME/nature_aml_log2_fpkm_cleaned.txt', sep='\t');head(exp_mat)

## drug sen data
drug_sen_df = as.data.frame(fread('../../data/VIZOME/nature_aml_drug_sen.tsv',sep='\t'))
drug_sen_df$lab_id = sapply(drug_sen_df$lab_id, function(x) str_replace_all(paste('X',x,sep=''),'-','.'))
drug_sen_df$lab_id = as.factor(as.character(drug_sen_df$lab_id))
drug_ls = drug_sen_df[['inhibitor']][!duplicated(drug_sen_df[['inhibitor']])]

## mut
NPM1_sample_ls = (maf_df%>%filter(Hugo_Symbol=='NPM1'))$Tumor_Sample_Barcode
NPM1_WT_sample_ls = setdiff(maf_df$Tumor_Sample_Barcode,NPM1_sample_ls)
NPM1_sample_ls = str_replace_all(paste('X',NPM1_sample_ls,sep=''),'-','.')
NPM1_WT_sample_ls = str_replace_all(paste('X',NPM1_WT_sample_ls,sep=''),'-','.')

NPM1_df = data.frame( sample=c(NPM1_sample_ls,NPM1_WT_sample_ls),
                     NPM1_status=as.factor(c(rep('mut',length(NPM1_sample_ls)), 
                                             rep('wt',length(NPM1_WT_sample_ls)))) )

# clinc
cli_df = as.data.frame(fread('/mnt/d/my_project/data/VIZOME/nature_aml_clinical_annotation.txt'))

blastbm_df = cli_df[c('LabId','%.Blasts.in.BM')]
blastbm_df$LabId = str_replace_all(paste('X',blastbm_df$LabId,sep=''),'-','.')
colnames(blastbm_df) = c('lab_id','Blasts.in.BM')
blastbm_df$Blasts.in.BM = as.numeric(blastbm_df$Blasts.in.BM)


our_drugs = c('Midostaurin','Quizartinib (AC220)','Crenolanib','Gilteritinib (ASP-2215)')

############
############ 
hl_8gene_set = strsplit("CEACAM8
    OLR1
    INHBA
    ANXA3
    FCAR
    S100A8
    PADI4
    SYNE1",'\n')[[1]]

raw_cebpa_set = strsplit("PRTN3	PRRG4	TRIB1	SLC51A	MEFV	SIGLEC16	CD1D	DHRS9	KCNE1	ITGB2	GRAMD1B	NCF2	GPR141	S100A8	S100A9	SNRPN	GPR160	SRGN	CCR2	AP5B1	CLEC5A	IPCEF1	ATG16L2	UBXN2B	SIGLEC9	ATG7	PNOC	GPR183	DENND1A	MFSD2A	IPMK	BCL2A1	EVI2B	NINL	STX10	FAM107B	UTRN	CAMP	KIAA0319L	ANXA4	ATP6V1E2	PTK2B	ST14	ADAP1	ALAS1	PIWIL4	APBB1IP	TMEM150B	FAM200B	UBAC2	ZNF438	SIRPB2	AIM2	HNRNPF	FPR2	ECHDC1	FPR1	RNLS	GBA	CCL3	FIG4	SLC35A1	TPD52	IER2	RCBTB2	BLM	FBXL5	BEST1	STK17A	LAIR1	ALOX5AP	PLD1	USP3	FAR2	ELANE	XPO6	ERLIN1	DAGLB	CISD1	POLD1	TRPM2	FOS	S100A11	IL10RB	NR1H2	IL6R	AZU1	NLRP3	SLTM	MARCO	PLEK	SMAP2	SERPINB8	IL1RAP	MKNK1	MR1	ITGAL	RCAN1	HDAC9	PTPN12	ACSL1	TMCO6	ANKRD11	PRDM10	QSOX1	MBNL1	PWWP2A	NANS	IFNGR2	FMO4	GLIPR1	FTH1	CES1	CNOT7	SSH2	COQ10B	CSRNP2	IL4R	MAP3K8	TMPPE	PTPN2	TMBIM4	LGALS3	HIPK3	NAPSA	P2RY2	PRR5L	DBN1	SNTB1	HNRNPU	FLT3	VPS37A	NUDT3	GLRX	SMYD3	COPZ1	CDC123	HIP1	PTPRE	SPRY4	MCOLN2	CDKN1A	SCARB1	ARNT	MTHFD2	DRAM1	NFE2	ST3GAL6	CTSC	RHOBTB2	NOD2	SRD5A3	SLC22A16	HIPK2	ANXA1	TNPO1	IKBKE	SLC25A37	SLC24A4	USP32	HLX	CASS4	ABCF2	ERI1	ZFP36L2	RTN4	GMEB2	TMEM104	IRF2BP2	GLB1	ST3GAL1	SPRED2	PIK3CB	NR3C1	TRPC4AP	ZC3H7A	ATG13	B3GNT2	RHOQ	STX6	PAQR8	TGFBR2	SPATA18	SMARCA2	NUDT5	ATP6V0A1	SLC17A5	PRKAR2B	ZNF395	NAT9	TPT1	PMAIP1	ARL5C	RAPGEF2	FBXO30	C12orf10	MED13L	ELMSAN1	CNST	LPP	MILR1	SEPHS2	PAFAH2	NREP	UNC79	AK9	LPAR1	EHD4	ZDHHC14	STARD13	HOXA5",
                        '\t')[[1]]

if(TRUE){new_cebpa_set = strsplit("GPR141
    GPR160
    BCL2A1
    S100A8
    AP5B1
    PAFAH2
    CD1D
    ITGB2
    DHRS9
    ALOX5AP
    S100A9
    CASS4
    MEFV
    SIGLEC9
    IPCEF1
    CLEC5A
    DENND1A
    USP3
    KCNE1
    IFNGR2
    CCR2
    SIGLEC16
    CAMP",'\n')[[1]]}
############
############

Mido_corr_gene = cal_corr(exp_mat,our_drugs[1],drug_sen_df,raw_cebpa_set)$Mido_corr_gene
Quiz_corr_gene = cal_corr(exp_mat,our_drugs[2],drug_sen_df,raw_cebpa_set)$Mido_corr_gene
Quiz_corr_gene = cal_corr(exp_mat,our_drugs[3],drug_sen_df,raw_cebpa_set)$Mido_corr_gene
Quiz_corr_gene = cal_corr(exp_mat,our_drugs[4],drug_sen_df,raw_cebpa_set)$Mido_corr_gene




cebpa_set = unique(c(Mido_corr_gene[1:5],
                     Quiz_corr_gene[1:5],
                     Creno_corr_gene[1:5],
                     Gilter_corr_gene[1:5], 
                     new_cebpa_set))

rownames(Quiz_corr_df) = NULL
rownames(Creno_corr_df) = NULL
rownames(Gilter_corr_df) = NULL
rownames(Mido_corr_df) = NULL
cebpa_set = cebpa_set[order((Quiz_corr_df%>%filter(gene%in%cebpa_set))$p+
(Creno_corr_df%>%filter(gene%in%cebpa_set))$p+
(Gilter_corr_df%>%filter(gene%in%cebpa_set))$p+
(Mido_corr_df%>%filter(gene%in%cebpa_set))$p,decreasing=TRUE)]


cebpa_new_ls_df = data.frame(new_cebpa_set=new_cebpa_set,quiz=rep(0,length(new_cebpa_set)),gilter=rep(0,length(new_cebpa_set)),
                             creno=rep(0,length(new_cebpa_set)),mido=rep(0,length(new_cebpa_set)))

cebpa_new_ls_df[cebpa_new_ls_df$new_cebpa_set%in%Mido_corr_df[1:10,]$gene,]$mido = 1
cebpa_new_ls_df[cebpa_new_ls_df$new_cebpa_set%in%Quiz_corr_df[1:10,]$gene,]$quiz = 1
cebpa_new_ls_df[cebpa_new_ls_df$new_cebpa_set%in%Creno_corr_df[1:10,]$gene,]$creno = 1
cebpa_new_ls_df[cebpa_new_ls_df$new_cebpa_set%in%Gilter_corr_df[1:10,]$gene,]$gilter = 1


good_gene_ls = cebpa_new_ls_df[apply(cebpa_new_ls_df[,c(2:5)],1,sum)!=0,]$new_cebpa_set

# good_gene_ls = cebpa_set
good_auc_ls = c(0,0,0, 0.6666667)
precent = 0.3
# 由于gilter不好，针对gilter进行优化
for(jj in 1:nn){
    message(sprintf('%s >\n',jj))
    
    last_biggest = 0
    last_biggest_i = 0
    tmp_Gilter_corr_gene = setdiff(Gilter_corr_gene, good_gene_ls)
    for(i in 1:length(tmp_Gilter_corr_gene)){
        #if(i%%20==0){
        #    message('.')
        #}
        tmp_gene_ls = unique(c(good_gene_ls,tmp_Gilter_corr_gene[i])) # del one gene
        # if()

        tmp_auc_ls = c(0,0,0,0)
        
        #for(ii in 1:4){
        drug = our_drugs[4]
        #print(drug)
        ## choose sp sample
        itd_exp_drug_sample = intersect(intersect(ITD_sample_ls,colnames(exp_mat)), (drug_sen_df%>%filter(inhibitor==drug))$lab_id)
        get_sen_unsen_res = get_sen_unsen(drug_sen_df=drug_sen_df, drug=drug,
                            sp_sample_ls=itd_exp_drug_sample, exp_mat=exp_mat, auc_ic50='auc',precent=precent)
        res_df_cebpa = test_gene_ls_in_drug_for_Exhaustion(drug_sen_df=drug_sen_df, drug=drug,
                            sp_sample_ls=itd_exp_drug_sample, sen_unsen_list=get_sen_unsen_res, 
                            exp_mat=exp_mat, gene_ls=tmp_gene_ls,auc_ic50='auc')
        roc1_2 <- roc(res_df_cebpa$df_plot$sen_unsen, 
                        res_df_cebpa$df_plot$ssGSEA_score,
                        quiet=TRUE)
        tmp_auc_ls[4] = roc1_2$auc
        if(last_biggest<roc1_2$auc){
            last_biggest = roc1_2$auc
            last_biggest_i = i
        }    

        #}
        #message(ii,' ',good_gene_ls[i],' ',paste(tmp_auc_ls,collapse=' '))
        if(good_auc_ls[4]<tmp_auc_ls[4]){

            for(ii in 1:3){
                drug = our_drugs[ii]
                #print(drug)
                ## choose sp sample
                itd_exp_drug_sample = intersect(intersect(ITD_sample_ls,colnames(exp_mat)), (drug_sen_df%>%filter(inhibitor==drug))$lab_id)
                get_sen_unsen_res = get_sen_unsen(drug_sen_df=drug_sen_df, drug=drug,
                                    sp_sample_ls=itd_exp_drug_sample, exp_mat=exp_mat, auc_ic50='auc',precent=precent)
                res_df_cebpa = test_gene_ls_in_drug_for_Exhaustion(drug_sen_df=drug_sen_df, drug=drug,
                                    sp_sample_ls=itd_exp_drug_sample, sen_unsen_list=get_sen_unsen_res, 
                                    exp_mat=exp_mat, gene_ls=tmp_gene_ls,auc_ic50='auc')
                roc1_2 <- roc(res_df_cebpa$df_plot$sen_unsen, 
                                res_df_cebpa$df_plot$ssGSEA_score,
                                quiet=TRUE)
                tmp_auc_ls[ii] = roc1_2$auc
            }
            message(i,' ',tmp_Gilter_corr_gene[i],' ',paste(tmp_auc_ls,collapse=' '))
            if(0.6900585<=tmp_auc_ls[1]&0.75<tmp_auc_ls[2]&0.75<tmp_auc_ls[3]){
                good_auc_ls[4] = tmp_auc_ls[4]
                good_gene_ls = tmp_gene_ls
                cat(good_gene_ls,sep=',')
                break                
            }
        }else{
            message(i,' ',tmp_Gilter_corr_gene[i],' ',paste(tmp_auc_ls,collapse=' '))
        }
            

    }
    if(0.75<good_auc_ls[4]){
        break
    }
    if(last_biggest<=good_auc_ls[4]){
        good_gene_ls = c(good_gene_ls, tmp_Gilter_corr_gene[last_biggest_i])
        Gilter_corr_gene = setdiff(Gilter_corr_gene,tmp_Gilter_corr_gene[last_biggest_i])
    }
}










