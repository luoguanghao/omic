
# 每个药物组合给一个热图展示
if(FALSE){

    path_list = list(pho=list(path='./data/drug_pho_th0.4_cor.xlsx',strong_positive_th=0.4),
                    pro=list(path='./data/drug_pro_th0.5_cor.xlsx',strong_positive_th=0.5),
                    rna=list(path='./data/drug_rna_th0.4_cor.xlsx',strong_positive_th=0.4))

    result_dir = sprintf('drug_ranking_20211122')
    dir.create(result_dir)

    drug_pair_df_list = list()
    for(term in names(path_list)){
        # term = 'pho'
        drug_comb_df = readxl::read_excel(path_list[[term]][['path']])
        colnames(drug_comb_df)[1] = 'feature'
        
        strong_positive_th = path_list[[term]][['strong_positive_th']]

        comb_res = combn(colnames(drug_comb_df)[-1],2)

        long_df = data.frame()
        drug_pair_df = list(comb_no=c(),pair=c(),dA=c(),dB=c(),pos=c(),neg=c(),score=c())
        drug_pair_ls = c()
        for(i in 1:dim(comb_res)[2]){
        #for(i in 1:5){
            drug_pair = comb_res[,i]
            tmp_df = drug_comb_df[,c('feature',drug_pair)]
            colnames(tmp_df)[2:3] = c('d1','d2')

            ttmp_df = tmp_df%>%filter((d1>.0 & d2<(-.0)) | (d2>.0&d1<(-.0))) # inverse correlation
            ttmp_df = ttmp_df%>%filter(d1>strong_positive_th | d2>strong_positive_th) # at least one strong positive correlation    

            if(dim(ttmp_df)[1]>=1){
                drug_pair_df$comb_no = c(drug_pair_df$comb_no, i)
                drug_pair_df$pair = c(drug_pair_df$pair, paste(drug_pair,collapse=' & '))
                drug_pair_df$dA = c(drug_pair_df$dA, drug_pair[1])
                drug_pair_df$dB = c(drug_pair_df$dB, drug_pair[2])
                drug_pair_df$pos = c(drug_pair_df$pos, dim((ttmp_df%>%filter(d1>0)))[1])
                drug_pair_df$neg = c(drug_pair_df$neg, dim((ttmp_df%>%filter(d1<0)))[1])
                score = (dim((ttmp_df%>%filter(d1<0)))[1]) * (dim((ttmp_df%>%filter(d1>0)))[1]) # A*B
                drug_pair_df$score = c(drug_pair_df$score, score)
            }
        }
        drug_pair_df = data.frame(drug_pair_df)
        drug_pair_df = drug_pair_df[order(drug_pair_df$score,decreasing=TRUE),]


        drug_pair_df_list[[term]] = drug_pair_df
        write_tsv(drug_pair_df, file=sprintf('%s/drug_ranking_%s_th%s.tsv',result_dir,term,strong_positive_th))
    }

}



if(FALSE){ # 早期朴素代码

    # example notebook: ~/my_project/hl/大组学药物联用/drug_ranking.ipynb
    library(tidyverse)


    LZW_drug_comb_df = readxl::read_excel('./data/LZW_drug_comb2.xlsx')
    colnames(LZW_drug_comb_df)[1] = 'phosph_site'
    # LZW_drug_comb_df = column_to_rownames(LZW_drug_comb_df,'phosph_site')





    comb_res = combn(colnames(LZW_drug_comb_df)[-1],2)

    long_df = data.frame()
    drug_pair_df = list(comb_no=c(),pair=c(),dA=c(),dB=c(),pos=c(),neg=c(),score=c())
    drug_pair_ls = c()
    for(i in 1:dim(comb_res)[2]){
    #for(i in 1:5){
        drug_pair = comb_res[,i]
        tmp_df = LZW_drug_comb_df[,c('phosph_site',drug_pair)]
        colnames(tmp_df)[2:3] = c('d1','d2')
        
        ttmp_df = tmp_df%>%filter((d1>.0 & d2<(-.0)) | (d2>.0&d1<(-.0))) # inverse correlation
        ttmp_df = ttmp_df%>%filter(d1>.5 | d2>.5) # at least one strong positive correlation    

        if(dim(ttmp_df)[1]>=1){
            drug_pair_df$comb_no = c(drug_pair_df$comb_no, i)
            drug_pair_df$pair = c(drug_pair_df$pair, paste(drug_pair,collapse=' & '))
            drug_pair_df$dA = c(drug_pair_df$dA, drug_pair[1])
            drug_pair_df$dB = c(drug_pair_df$dB, drug_pair[2])
            drug_pair_df$pos = c(drug_pair_df$pos, dim((ttmp_df%>%filter(d1>0)))[1])
            drug_pair_df$neg = c(drug_pair_df$neg, dim((ttmp_df%>%filter(d1<0)))[1])
            score = (dim((ttmp_df%>%filter(d1<0)))[1]) * (dim((ttmp_df%>%filter(d1>0)))[1]) # A*B
            drug_pair_df$score = c(drug_pair_df$score, score)
        }
    }
    drug_pair_df = data.frame(drug_pair_df)
    drug_pair_df = drug_pair_df[order(drug_pair_df$score,decreasing=TRUE),]



    write_tsv(drug_pair_df, file='./drug_ranking2.tsv')


}




