

#############
# rnaseq
#############
exp_dir = '/y/Archive/Bondi/data/VIZOME/nature_aml_log2_fpkm.txt'
###
probeMap_df = read.csv(probeMap_dir,sep='\t')
exp_mat = read.csv(exp_dir,sep='\t')

exp_mat = aggregate(.~ Symbol, exp_mat, mean)
exp_mat = exp_mat %>% column_to_rownames('Symbol')


### output
# exp_mat
# sample_ls = colnames(exp_mat)
#############
#############



#############
# mut
#############
maf_dir = '/y/Archive/Bondi/data/VIZOME/LAML_nature_PDC.maf'
###
maf_df = read.csv(maf_dir, sep='\t')

### 获取FLT3_ITD
ITD_sample_ls = as.vector((maf_df %>% filter(Hugo_Symbol=='FLT3'&ITDorNOT=='ITD'))['Tumor_Sample_Barcode'][['Tumor_Sample_Barcode']])
ITD_sample_ls = str_replace_all(paste('X',ITD_sample_ls,sep=''),'-','.') # colname的格式化
###
# maf_df
# ITD_sample_ls
#############
#############




#############
# drug sen
#############
drugsen_dir = '/y/Archive/Bondi/data/VIZOME/nature_aml_drug_sen.tsv'
###
drug_sen_df = read.csv(drugsen_dir, sep='\t')

drug_sen_df$lab_id = lapply(drug_sen_df$lab_id, function(x) str_replace_all(paste('X',x,sep=''),'-','.')) # colname的格式化
drug_sen_df$lab_id = as.factor(as.character(drug_sen_df$lab_id))
drug_ls = drug_sen_df[['inhibitor']][!duplicated(drug_sen_df[['inhibitor']])] 
###
# drug_sen_df # every drug every sample, a long long frame
# drug_ls
#############
#############



#############
# sample info
#############
info_dir = '/y/Archive/Bondi/data/VIZOME/nature_aml_drug_sen.tsv'
###
info_df = read.csv(info_dir, sep='\t')


###
# drug_sen_df # every drug every sample, a long long frame
# drug_ls
#############
#############












