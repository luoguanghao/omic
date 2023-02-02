library(tidyverse)

######## ch pro
phogene_mat = read_tsv('/home/lgh/my_project/data/changhai_data/protein_phosph/second/sort_by_sample/impute/new/phosphogene_df_min_impute_all.tsv')%>%column_to_rownames('Hugo_Symbol')
pho_mat = read_tsv('/home/lgh/my_project/data/changhai_data/protein_phosph/second/sort_by_sample/impute/new/phosphosite_df_min_impute_all.tsv')%>%column_to_rownames('Hugo_Symbol')
pro_mat = read_tsv('/home/lgh/my_project/data/changhai_data/protein_phosph/second/sort_by_sample/impute/new/protein_df_min_impute_all.tsv')%>%column_to_rownames('Hugo_Symbol')

pro_mat_t = scale(t(pro_mat))

phogene_mat_raw = read_tsv('/home/lgh/my_project/data/changhai_data/protein_phosph/second/sort_by_sample/CH_amlid_phosphoprotein_gene_withna_nolog_geneXsample.tsv')%>%column_to_rownames('Hugo_Symbol')
pho_mat_raw = read_tsv('/home/lgh/my_project/data/changhai_data/protein_phosph/second/sort_by_sample/CH_amlid_phosphoprotein_site_withna_nolog_geneXsample.tsv')%>%column_to_rownames('Hugo_Symbol')
pro_mat_raw = read_tsv('/home/lgh/my_project/data/changhai_data/protein_phosph/second/sort_by_sample/CH_amlid_protein_withna_nolog_geneXsample.tsv')%>%column_to_rownames('Hugo_Symbol')

# sel_progene_by_na = rownames(phogene_mat_raw)[apply(phogene_mat_raw,1,function(x) sum(is.na(x)))<101*0.5]
sel_pro_by_na = rownames(pro_mat_raw)[apply(pro_mat_raw,1,function(x) sum(is.na(x)))<101*0.5]
sel_pho_by_na = rownames(pho_mat_raw)[apply(pho_mat_raw,1,function(x) sum(is.na(x)))<101*0.5]
sel_progene_by_na = unique(sapply(strsplit((sel_pho_by_na),'_'),'[[',1))

ch_clinical_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/clinical.tsv')


######
ch_exp_count_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_deseq2norm_expression_cleaned.tsv')
ch_exp_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_deseq2norm_expression_cleaned.tsv')
ch_exp_tpm_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_sampleid_tpm_expression_cleaned.tsv')
ch_surv_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_surv_rfs_eln_fab_0401_cleaned.tsv')
ch_exp_mat_scale = t(scale(t(ch_exp_df%>%column_to_rownames('Hugo_Symbol'))))

# ch_pro_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_tpm_expression_cleaned.tsv')
ch_drug_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_rawic50_drug_sensitivity_cleaned.tsv')
ch_mut_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_mutations_matrix_cleaned.tsv')


vizome_drug_anno_df = read_tsv('/home/lgh/my_project/data/VIZOME/cleaned/VIZOME_drug_anno_df.tsv')
vizome_drug_df = read_tsv('/home/lgh/my_project/data/VIZOME/cleaned/VIZOME_drug_sen_cleaned.tsv')
vizome_exp_df = read_tsv('/home/lgh/my_project/data/VIZOME/cleaned/VIZOME_log2fpkm_expression_cleaned.tsv')
vizome_surv_df = read_tsv('/home/lgh/my_project/data/VIZOME/cleaned/VIZOME_survival_cleaned.tsv')
vizome_exp_mat_scale = t(scale(t(vizome_exp_df%>%column_to_rownames('Hugo_Symbol'))))
vizome_clinic_df = read_tsv('/home/lgh/my_project/data/VIZOME/cleaned/clinical.tsv')


###### TCGA
exp_mat_tcga = read_tsv('/home/lgh/my_project/data/TCGA/LAML/TCGA-LAML.htseq_fpkm_cleaned.tsv')%>%column_to_rownames('Hugo_Symbol')
exp_mat_tcga_scale = t(na.omit(t(scale(t(exp_mat_tcga)))))
surv_tcga_df = read_tsv('/home/lgh/my_project/data/TCGA/LAML/TCGA-LAML.survival.cleaned.tsv')
surv_tcga_df$sample = str_replace_all(surv_tcga_df$sample,'[^[:alnum:]]','.')
tcga_clinic_df = read_tsv('/home/lgh/my_project/data/TCGA/LAML/TCGA-LAML.GDC_phenotype.tsv.gz')

tcga_itd_df = read_tsv('/home/lgh/my_project/data/TCGA/LAML/tcga_itd.tsv')
tcga_itd_df$sample = sapply(tcga_itd_df$sample, function(x) substr(x,1,(nchar(x)-4)))


######## AML2022 
pro_2022_mat = read_csv('/home/lgh/my_project/data/AML_Proteomics_2022/pro_impute_Discovery_Cohort.csv')
pro_2022_mat = cbind(pro_2022_mat[-1,179],pro_2022_mat[-1,-c(178:181)])

pro_2022_mat[,2:ncol(pro_2022_mat)] = sapply(pro_2022_mat[,2:ncol(pro_2022_mat)],as.numeric)
pro_2022_mat = aggregate(.~PG.Genes,pro_2022_mat,mean)
pro_2022_mat[[1]] = sapply(pro_2022_mat[[1]], function(x) if(length(strsplit(x,'\\.')[[1]])>1){strsplit(x,'\\.')[[1]][2]}else{x} )
colnames(pro_2022_mat)[1] = 'Hugo_Symbol'
pro_2022_mat_t = data.frame(t(pro_2022_mat%>%column_to_rownames('Hugo_Symbol')))
# colnames(pro_2022_mat_t) = sapply(colnames(pro_2022_mat_t), function(x) if(length(strsplit(x,'\\.')[[1]])>1){strsplit(x,'\\.')[[1]][2]}else{x} )
pro_2022_mat_t = scale(pro_2022_mat_t)
                                  
pro2022_survival = read_tsv('/home/lgh/my_project/data/AML_Proteomics_2022/Clinical_Discovery_Cohort.txt')[,c('ID','Death Event','OS [months]','Relapse Event','RFS [months]')]
colnames(pro2022_survival) = c('sample','OSS','OS','RFSS','RFS')


pro_2022_valid_mat = read_csv('/home/lgh/my_project/data/AML_Proteomics_2022/pro_impute_Validation_Cohort.csv',skip=1)
pro_2022_valid_mat = cbind(pro_2022_valid_mat[-1,91],pro_2022_valid_mat[-1,-c(76:91)])
colnames(pro_2022_valid_mat)[1] = 'Gene names'
pro_2022_valid_mat[,2:ncol(pro_2022_valid_mat)] = sapply(pro_2022_valid_mat[,2:ncol(pro_2022_valid_mat)],as.numeric)
pro_2022_valid_mat = aggregate(.~`Gene names`,pro_2022_valid_mat,mean)
pro_2022_valid_mat[['Gene names']] = sapply(pro_2022_valid_mat[['Gene names']], function(x) if(length(strsplit(x,'\\.')[[1]])>1){strsplit(x,'\\.')[[1]][2]}else{x} )
colnames(pro_2022_valid_mat)[1] = 'Hugo_Symbol'
pro_2022_valid_mat_t = data.frame(t(pro_2022_valid_mat%>%column_to_rownames('Hugo_Symbol')))
# colnames(pro_2022_valid_mat_t) = sapply(colnames(pro_2022_valid_mat_t), function(x) if(length(strsplit(x,'\\.')[[1]])>1){strsplit(x,'\\.')[[1]][2]}else{x} )

pro2022_valid_survival = read_tsv('/home/lgh/my_project/data/AML_Proteomics_2022/Clinical_Validation_Cohort.txt')[,c('ID','Death Event','OS [months]','Relapse Event','RFS [months]')]
colnames(pro2022_valid_survival) = c('sample','OSS','OS','RFSS','RFS')
                                        
pro_2022_valid_mat_t = scale(pro_2022_valid_mat_t)


pro_2022_mut = read_tsv('/home/lgh/my_project/data/AML_Proteomics_2022/Clinical_Discovery_Cohort.txt')[,c('ID','NPM1','FLT3')]
pro_2022_valid_mut = read_tsv('/home/lgh/my_project/data/AML_Proteomics_2022/Clinical_Validation_Cohort.txt')[,c('ID','NPM1','FLT3')]
colnames(pro_2022_mut)[1] = 'sample'
colnames(pro_2022_valid_mut)[1] = 'sample'


############## gmt

PATH_DATA = '/home/lgh/my_project/data'
gmtfile = list(
    gobp=paste(PATH_DATA,'/pathway/c5.go.bp.v7.5.1.symbols.gmt',sep='')    
    #kegg=paste(PATH_DATA,'/pathway/c2.cp.kegg.v7.5.1.symbols.gmt',sep='')
    #reactome=paste(PATH_DATA,'/pathway/c2.cp.reactome.v7.5.1.symbols.gmt',sep=''),
    #hallmark='/home/lgh/my_project/data/pathway/h.all.v7.5.1.symbols.gmt'
    #tft=paste(PATH_DATA,'/pathway/c3.tft.v7.5.1.symbols.add_cebpa.gmt',sep='')
)













