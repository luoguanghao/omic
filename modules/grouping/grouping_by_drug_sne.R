

# get_sen_unsen_sample
## sp_drug_sen_df: in vizome format:inhibitor;lab_id;auc;ic50
## example: sp_drug_sen_df = drug_sen_df[drug_sen_df['inhibitor']==drug,]%>%filter(lab_id%in%sp_sample_ls)
get_sen_unsen_sample <- function(sp_drug_sen_df, by='auc', precent=0.3,){

    sp_drug_sen_df = sp_drug_sen_df[order(sp_drug_sen_df$auc),]
    sp_drug_sen_df$lab_id = factor(sp_drug_sen_df$lab_id,level=as.vector(sp_drug_sen_df$lab_id))

    top_bottom_number = round(dim(sp_drug_sen_df)[1]*0.3)
    sen_sample = as.vector(sp_drug_sen_df[1:top_bottom_number,]$lab_id)
    unsen_sample = as.vector(sp_drug_sen_df[(dim(sp_drug_sen_df)[1]-top_bottom_number):dim(sp_drug_sen_df)[1],]$lab_id)

    return(list(sen_sample_ls = sen_sample, unsen_sample_ls = unsen_sample))
}


















