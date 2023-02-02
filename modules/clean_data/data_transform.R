# https://cloud.tencent.com/developer/article/1526683


# wind2long data #

# 就是每列数据拿来拼接变成一长列数据
# tidy::gather
long_data = data %>% gather(key = "variables", value = "values")

# reshape2::melt
## id.vars:不需要melt的列
#
# melt(data, id.vars, measure.vars,
#   variable.name = "variable", ..., na.rm = FALSE,
#   value.name = "value")
#
data3 <- reshape2::melt(iris, id.vars = c("Species"))

# long2wind #
## spread or dcast
reshape2::dcast(cal_mean, inhibitor~group, value.var='median'))


drug_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_log10ic50_drug_sensitivity_cleaned.tsv')
reshape2::dcast(drug_df, inhibitor~sample, value.var='ic50')






#########################
#########################
mut_df = read_tsv('/home/lgh/my_project/data/changhai_data/cleaned/CH_amlid_mutations_cleaned.maf')

mut_mat_df = data.frame(matrix(rep(0,length(unique(mut_df$sample))*length(unique(mut_df$Hugo_Symbol))),
      nrow = length(unique(mut_df$sample)), ncol = length(unique(mut_df$Hugo_Symbol))))
rownames(mut_mat_df) = unique(mut_df$Tumor_Sample_Barcode)
colnames(mut_mat_df) = unique(mut_df$Hugo_Symbol)
mut_mat_df = mut_mat_df%>%rownames_to_column('sample')

for(i in 2:ncol(mut_mat_df)){
    mut_sample_vec = (mut_df%>%filter(Hugo_Symbol==colnames(mut_mat_df)[i]))$Tumor_Sample_Barcode
    mut_mat_df[mut_mat_df$sample%in%mut_sample_vec,][[i]] = 1
}
mut_mat_df[,c(2:ncol(mut_mat_df))] = sapply(mut_mat_df[,c(2:ncol(mut_mat_df))],as.factor)
# mut_mat_df

drug='ABT-263 (Navitoclax)'
df_for_cal = merge((drug_sen_df%>%filter(inhibitor==drug))[,c('sample','inhibitor','ic50')],mut_mat_df,by='sample')














