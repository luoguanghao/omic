#############
# TCGA Data
#############
'''
TCGA format
exp_mat = gene_id * sample
probeMap: gene_id;gene_symbol;other_info
'''

exp_dir = '~/my_project/data/TCGA/LAML/TCGA-LAML.htseq_fpkm.tsv.gz'
probeMap_dir = '~/my_project/data/TCGA/LAML/gencode.v22.annotation.gene.probeMap'
######

exp_mat = pd.read_csv(exp_dir,sep='\t')
probeMap_df = pd.read_csv(probeMap_dir,sep='\t')

# do anno;format the mat;remove dup gene
exp_mat = pd.merge(probeMap_df.iloc[:,0:2],exp_mat,left_on='id',right_on='Ensembl_ID').drop(['id','Ensembl_ID'],axis=1)

######
exp_mat.loc[exp_mat.duplicated(keep=False)==False].to_csv('LAML_rpkm_symbol_nodup.tsv',sep='\t')
# return exp_mat
#############
#############


#############
# Vizome Data
#############



#############
#############






#############
# Changhai Data
#############




#############
#############









