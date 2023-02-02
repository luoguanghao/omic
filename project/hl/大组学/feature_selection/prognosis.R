
library(tidyverse)
library(fastcox)
library(timeROC)
source('/home/lgh/my_project/omic/modules/plot/survival.R')





pho_mat = read_tsv('/home/lgh/my_project/data/changhai_data/protein_phosph/second/sort_by_sample/impute/new/phosphosite_df_min_impute_scalegene.tsv')%>%column_to_rownames('Hugo_Symbol')
pro_mat = read_tsv('/home/lgh/my_project/data/changhai_data/protein_phosph/second/sort_by_sample/impute/new/protein_df_min_impute_scalegene.tsv')%>%column_to_rownames('Hugo_Symbol')

pho_mat_raw = read_tsv('/home/lgh/my_project/data/changhai_data/protein_phosph/second/sort_by_sample/CH_amlid_phosphoprotein_site_withna_nolog_geneXsample.tsv')%>%column_to_rownames('Hugo_Symbol')
pro_mat_raw = read_tsv('/home/lgh/my_project/data/changhai_data/protein_phosph/second/sort_by_sample/CH_amlid_protein_withna_nolog_geneXsample.tsv')%>%column_to_rownames('Hugo_Symbol')

sel_pro_by_na = rownames(pro_mat_raw)[apply(pro_mat_raw,1,function(x) sum(is.na(x)))<101*0.7]
sel_pho_by_na = rownames(pho_mat_raw)[apply(pho_mat_raw,1,function(x) sum(is.na(x)))<101*0.7]  

useful_sample = intersect(colnames(pro_mat),(surv_df%>%filter(!is.na(OS)))$sample)



pro_2022_mat = read_csv('/home/lgh/my_project/data/AML_Proteomics_2022/pro_impute_Discovery_Cohort.csv')
pro_2022_mat = cbind(pro_2022_mat[-1,179],pro_2022_mat[-1,-c(178:181)])

pro_2022_mat[,2:ncol(pro_2022_mat)] = sapply(pro_2022_mat[,2:ncol(pro_2022_mat)],as.numeric)
pro_2022_mat = aggregate(.~PG.Genes,pro_2022_mat,mean)
pro_2022_mat_t = data.frame(t(pro_2022_mat%>%column_to_rownames('PG.Genes')))
colnames(pro_2022_mat_t) = sapply(colnames(pro_2022_mat_t), function(x) if(length(strsplit(x,'\\.')[[1]])>1){strsplit(x,'\\.')[[1]][2]}else{x} )


pro2022_survival = read_tsv('/home/lgh/my_project/data/AML_Proteomics_2022/Clinical_Discovery_Cohort.txt')[,c('ID','Death Event','OS [months]')]
colnames(pro2022_survival) = c('sample','OSS','OS')



# ===================
x = t(pro_mat[sel_pro_by_na,useful_sample])
y = (surv_df%>%column_to_rownames('sample'))[useful_sample,]$OS
status = (surv_df%>%column_to_rownames('sample'))[useful_sample,]$OSS

useful_gene_comb_pro2022 = intersect(colnames(pro_2022_mat_t),colnames(x))


# 迭代参数与seed，训练模型，测试集验证

res_list_comb_pro2022_no_divi = list()

for(seed_tmp in 1:25){
    message(sprintf('# %s',seed_tmp))
    tmp_res = list()

    set.seed(seed_tmp)
    sample_res = useful_sample[1:length(useful_sample)]

    #x_valid = x[sample_res, useful_gene_comb_pro2022]
    #y_valid = y[sample_res]
    #status_valid = status[sample_res]

    x_do = x[setdiff(1:length(useful_sample),sample_res), useful_gene_comb_pro2022]
    y_do = y[setdiff(1:length(useful_sample),sample_res)]
    status_do = status[setdiff(1:length(useful_sample),sample_res)]

    cv1<-cv.cocktail(x=x_do,y=y_do,d=status_do,alpha=0.5,nfolds=10) 

    lbd_seq = seq(cv1$lambda.min, cv1$lambda.1se, (cv1$lambda.1se-cv1$lambda.min)/10)
    
    if(length(lbd_seq)!=1){
    
        m1<-cocktail(x=x_do,y=y_do,d=status_do,alpha=0.5,lambda=lbd_seq)

        #selected_featrue = as.matrix(predict(m1,newx=x_do,type="coefficients"))[,'s40']
        #selected_featrue = selected_featrue[selected_featrue>0]

        res_df_1 = list(ldb=c(),p=c(),n_high=c(),n_low=c())
        df_cal_1 = merge(surv_df,data.frame(predict(m1,newx=x_do,type="response"))%>%rownames_to_column('sample'))
        surv_res_1 = list()
        for(i in 1:(length(lbd_seq)-1)) {
            surv_res = plot_survival_contin(df_cal_1, sprintf('s%s',i), plot=F, title=NULL, OS_RF='OS')

            res_df_1$ldb = c(res_df_1$ldb, sprintf('s%s',i))
            res_df_1$p = c(res_df_1$p, surv_res$pval)

            res_df_1$n_high = c(res_df_1$n_high, surv_res$surv_diff$n[[1]])
            res_df_1$n_low = c(res_df_1$n_low, surv_res$surv_diff$n[[2]])

            surv_res_1[[i]] = surv_res
            #print(plot_survival_contin(df_cal, sprintf('s%s',i), plot=T, title=NULL, OS_RF='OS')$plot+labs(title=sprintf('s%s',i)))
        }
        res_df_1 = data.frame(res_df_1)


        res_df_2 = list(ldb=c(),p=c(),n_high=c(),n_low=c())
        df_cal_2 = merge(pro2022_survival, 
                         data.frame(predict(m1,newx=as.matrix(scale(pro_2022_mat_t))[,useful_gene_comb_pro2022],type="response"))%>%
                             rownames_to_column('sample'))
        surv_res_2 = list()
        for(i in 1:(length(lbd_seq)-1) ){
            if(Inf %in% (df_cal_2[[sprintf('s%s',i)]])){
                res_df_2$ldb = c(res_df_2$ldb, sprintf('s%s',i))
                res_df_2$p = c(res_df_2$p, 'Inff')

                res_df_2$n_high = c(res_df_2$n_high, 'Inff')
                res_df_2$n_low = c(res_df_2$n_low, 'Inff')

                surv_res_2[[i]] = 'Inff'

            }else{
                surv_res = plot_survival_contin(df_cal_2, sprintf('s%s',i), plot=F, title=NULL, OS_RF='OS')
                #print(plot_survival_contin(df_cal, sprintf('s%s',i), plot=T, title=NULL, OS_RF='OS')$plot+labs(title=sprintf('s%s',i)))
                res_df_2$ldb = c(res_df_2$ldb, sprintf('s%s',i))
                res_df_2$p = c(res_df_2$p, surv_res$pval)

                res_df_2$n_high = c(res_df_2$n_high, surv_res$surv_diff$n[[1]])
                res_df_2$n_low = c(res_df_2$n_low, surv_res$surv_diff$n[[2]])

                surv_res_2[[i]] = surv_res
            }
        }
        res_df_2 = data.frame(res_df_2)

        #res_list_comb_pro2022[[ii]]$df_cal_pro2022 = df_cal_2
        #res_list_comb_pro2022[[ii]]$res_df_pro2022 = res_df_2
        #res_list_comb_pro2022[[ii]]$surv_res_pro2022 = surv_res_2



        tmp_res = list(
            sample_res=sample_res,
            cv1=cv1,
            m1=m1,
            #selected_featrue=selected_featrue,

            df_cal_1=df_cal_1,
            df_cal_pro2022=df_cal_2,

            res_df_1 = res_df_1,
            res_df_pro2022 = res_df_2,

            surv_res_1 = surv_res_1,
            surv_res_pro2022 = surv_res_2
        )
    }else{
        tmp_res = list(
            sample_res=sample_res,
            cv1=cv1,
            # m1=m1,
            #selected_featrue=selected_featrue,

            df_cal_1="NO",
            df_cal_pro2022="NO",

            res_df_1 = "NO",
            res_df_pro2022 = "NO",

            surv_res_1 = "NO",
            surv_res_pro2022 = "NO"
        )        
    }
    res_list_comb_pro2022_no_divi[[seed_tmp]] = tmp_res

}


# plot lasso plot

for(ii in 1:length(res_list_comb_pro2022_no_divi)){
    plot(res_list_comb_pro2022_no_divi[[ii]]$cv1,main=ii)
}


# plot lambda min lse for every seed
for(i in c(2, 5, 6, 7, 8, 10,11,12,14,15,19,20,21,23,24)){
    print(c(sprintf('%s',i),res_list_comb_pro2022_no_divi[[i]]$cv1$lambda.min, res_list_comb_pro2022_no_divi[[i]]$cv1$lambda.1se))
}


# 获得各个seed-lambda组合的预测结果 p值

tmp = res_list_comb_pro2022_no_divi[[2]]$res_df_pro2022
tmp$seed = 2
res_df_pro2022_no_divi_all = tmp
for(ii in 3:length(res_list_comb_pro2022_no_divi)){
    
    tmp = res_list_comb_pro2022_no_divi[[ii]]$res_df_pro2022
    if(tmp[[1]]=='NO') next
    
    tmp$seed = ii
    res_df_pro2022_no_divi_all = rbind(res_df_pro2022_no_divi_all, tmp)
    
}

###############
# KM检验与作图
###############

# cut median 中位数cut来作图
ttmp_df = res_df_pro2022_no_divi_all%>%filter(p<0.05)
p_tmp = 0
for(ii in 1:dim(ttmp_df)[1]){
    # res = plot_survival_contin(res_list_comb_pro2022_no_divi[[ttmp_df[ii,]$seed]]$df_cal_pro2022, ttmp_df[ii,]$ldb, plot=T, title=NULL, OS_RF='OS')
    tmp = res_list_comb_pro2022_no_divi[[ttmp_df[ii,]$seed]]$df_cal_1[,c('sample','OSS','OS',ttmp_df[ii,]$ldb)]
    tmp$group = 'high'
    tmp[tmp[[ttmp_df[ii,]$ldb]]<median(tmp[[ttmp_df[ii,]$ldb]]),]$group = 'low'
    print( plot_survival_binary(tmp, 'group', plot=T, conf.int=FALSE, OS_RF='OS')$plot+labs(title=sprintf('%s_%s',ttmp_df[ii,]$seed,ttmp_df[ii,]$ldb)) )
}

ttmp_df = res_df_pro2022_no_divi_all%>%filter(p<0.05)
p_tmp = 0
for(ii in 1:dim(ttmp_df)[1]){
    # res = plot_survival_contin(res_list_comb_pro2022_no_divi[[ttmp_df[ii,]$seed]]$df_cal_pro2022, ttmp_df[ii,]$ldb, plot=T, title=NULL, OS_RF='OS')
    tmp = res_list_comb_pro2022_no_divi[[ttmp_df[ii,]$seed]]$df_cal_pro2022[,c('sample','OSS','OS',ttmp_df[ii,]$ldb)]
    tmp$group = 'high'
    tmp[tmp[[ttmp_df[ii,]$ldb]]<median(tmp[[ttmp_df[ii,]$ldb]]),]$group = 'low'
    print( plot_survival_binary(tmp, 'group', plot=T, conf.int=FALSE, OS_RF='OS')$plot+labs(title=sprintf('%s_%s',ttmp_df[ii,]$seed,ttmp_df[ii,]$ldb)) )
}



# 最佳cutpoint作图
ttmp_df = res_df_pro2022_no_divi_all%>%filter(p<0.05)
p_tmp = 0
for(ii in 1:dim(ttmp_df)[1]){
    res = plot_survival_contin(res_list_comb_pro2022_no_divi[[ttmp_df[ii,]$seed]]$df_cal_1, ttmp_df[ii,]$ldb, plot=T, title=NULL, OS_RF='OS')
    if(p_tmp!=res$pval){
        print(res$plot+labs(title=sprintf('%s_%s',ttmp_df[ii,]$seed,ttmp_df[ii,]$ldb)))
        p_tmp=res$pval
    }
    # print(plot_survival_contin(res_list[[ttmp_df[ii,]$seed]]$df_cal_2, ttmp_df[ii,]$ldb, plot=T, title=NULL, OS_RF='OS')$plot)
}

ttmp_df = res_df_pro2022_no_divi_all
p_tmp = 0
for(ii in 1:dim(ttmp_df)[1]){
    res = plot_survival_contin(res_list_comb_pro2022_no_divi[[ttmp_df[ii,]$seed]]$df_cal_pro2022, ttmp_df[ii,]$ldb, plot=T, title=NULL, OS_RF='OS')
    #res = coxph(as.formula(sprintf('Surv(OS, OSS) ~ %s',ttmp_df[ii,]$ldb)), data = res_list_comb_pro2022_no_divi[[ttmp_df[ii,]$seed]]$df_cal_pro2022)
    print(res)
    print('###################')
    # print(plot_survival_contin(res_list[[ttmp_df[ii,]$seed]]$df_cal_2, ttmp_df[ii,]$ldb, plot=T, title=NULL, OS_RF='OS')$plot)
}



#################
# coxph检验与作图
#################

ttmp_df = res_df_pro2022_no_divi_all
p_tmp = 0
for(ii in 1:dim(ttmp_df)[1]){
    res = coxph(as.formula(sprintf('Surv(OS, OSS) ~ %s',ttmp_df[ii,]$ldb)), data = res_list_comb_pro2022_no_divi[[ttmp_df[ii,]$seed]]$df_cal_1)
    print(res)
    print('###################')
}

ttmp_df = res_df_pro2022_no_divi_all
p_tmp = 0
for(ii in 1:dim(ttmp_df)[1]){
    res = coxph(as.formula(sprintf('Surv(OS, OSS) ~ %s',ttmp_df[ii,]$ldb)), data = res_list_comb_pro2022_no_divi[[ttmp_df[ii,]$seed]]$df_cal_pro2022)
    print(res)
    print('###################')
}


#########################################
# 各个seed下 每个lambda下的feature，df情况
#########################################

for(ii in 1:25) {print(ii);print(res_list_comb_pro2022_no_divi[[ii]]$m1)}


#################
# 提取feature
#################

selected_featrue = as.matrix(predict(res_list_comb_pro2022_no_divi[[23]]$m1,newx=x_do,type="coefficients"))[,'s10']
selected_featrue = selected_featrue[abs(selected_featrue)>0]





#################
# plot heatmap
#################

options(repr.plot.width=10, repr.plot.height=15)
pheatmap::pheatmap(
    x_do[,names(selected_featrue)],
    cluster_cols = T,
    #cluster_rows = FALSE,
    scale = 'row',
    color = colorRampPalette(c("green","black", "red"))(50),
    # filename='./cencluster_hmap1.png',
    cellwidth = , cellheight = 10,border=FALSE,
    #cutree_cols=2, 
    #annotation_col  = anno_nc_df,
    #main = drug
)



#################
# plot line
#################

eln_pro = readxl::read_excel('eln蛋白.xlsx')

data3 <- reshape2::melt(eln_pro[,c('gene','median_ ELN1','median_ ELN2','median_ ELN3')], id.vars = c("gene"))
data3$variable = as.vector(data3$variable)
data3 = data3%>%filter(gene%in%names(selected_featrue))

data3[data3$variable=='median_ ELN1',]$variable='ELN1'
data3[data3$variable=='median_ ELN2',]$variable='ELN2'
data3[data3$variable=='median_ ELN3',]$variable='ELN3'
data3$variable = factor(data3$variable,level=paste('ELN',1:3,sep=''))
data3$value = as.numeric(data3$value)

options(repr.plot.width=6, repr.plot.height=5)
p = ggplot(data = data3)# + 
for(gene_tmp in unique(data3[[1]])){
    
    p = p+geom_line(data=data3%>%filter(gene==gene_tmp), mapping = aes(x = variable, y = value, group=1,color =gene))
    
}

p = p + geom_vline(xintercept=1:3, color='grey')
# data3$variable = as.factor(data3$variable)
p + scale_y_continuous(limits=c(0.5, 1.5), breaks=seq(0.5,1.5,0.5)) + 
    xlab('ELN')+ylab('Intensity')+
    #theme_classic() +
    theme(panel.grid.major =element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),#去除背景
                   panel.border = element_blank()) #去除边框





#################
# ROC
#################


mayoscore4 <- timeROC(T = res_list_comb_pro2022_no_divi[[23]]$df_cal_1$OS, 
                      delta = res_list_comb_pro2022_no_divi[[23]]$df_cal_1$OSS, 
                      marker = res_list_comb_pro2022_no_divi[[23]]$df_cal_1$s10,
    cause = 1, weighting = "marginal", times = c(365 * 1,365 * 2, 365 * 3, 365 * 4), 
                      ROC = TRUE, iid=TRUE)
plot(mayoscore4, time = 365 * 3, col = "#1c61b6", lwd = 2)
plotAUCcurve(mayoscore4, FP = 2, add = FALSE, conf.int = TRUE, conf.band = FALSE,
    col = "#1c61b6")




mayoscore4 <- timeROC(T = res_list_comb_pro2022_no_divi[[23]]$df_cal_pro2022$OS, 
                      delta = res_list_comb_pro2022_no_divi[[23]]$df_cal_pro2022$OSS, 
                      marker = res_list_comb_pro2022_no_divi[[23]]$df_cal_pro2022$s10,
    cause = 1, weighting = "marginal", times = c(12 * 1,12 * 2, 12 * 3, 12 * 4, 12 * 5, 12 * 6, 12 * 7), 
                      ROC = TRUE, iid=TRUE)
plot(mayoscore4, time = 12 * 4, col = "#1c61b6", lwd = 2)
plotAUCcurve(mayoscore4, FP = 2, add = FALSE, conf.int = TRUE, conf.band = FALSE,
    col = "#1c61b6")


###################
###################
# valid pro2022
###################
###################

pro_2022_valid_mat=read_csv('/home/lgh/my_project/data/AML_Proteomics_2022/pro_impute_Validation_Cohort.csv',skip=1)

pro_2022_valid_mat = read_csv('/home/lgh/my_project/data/AML_Proteomics_2022/pro_impute_Validation_Cohort.csv',skip=1)
pro_2022_valid_mat = cbind(pro_2022_valid_mat[-1,91],pro_2022_valid_mat[-1,-c(76:91)])
colnames(pro_2022_valid_mat)[1] = 'Gene names'
pro_2022_valid_mat[,2:ncol(pro_2022_valid_mat)] = sapply(pro_2022_valid_mat[,2:ncol(pro_2022_valid_mat)],as.numeric)
pro_2022_valid_mat = aggregate(.~`Gene names`,pro_2022_valid_mat,mean)
pro_2022_valid_mat_t = data.frame(t(pro_2022_valid_mat%>%column_to_rownames('Gene names')))
colnames(pro_2022_valid_mat_t) = sapply(colnames(pro_2022_valid_mat_t), function(x) if(length(strsplit(x,'\\.')[[1]])>1){strsplit(x,'\\.')[[1]][2]}else{x} )
pro_2022_valid_mat_t = scale(pro_2022_valid_mat_t)

pro2022_valid_survival = read_tsv('/home/lgh/my_project/data/AML_Proteomics_2022/Clinical_Validation_Cohort.txt')[,c('ID','Death Event','OS [months]')]
colnames(pro2022_valid_survival) = c('sample','OSS','OS')


add_mat = matrix(0,nrow=nrow(pro_2022_valid_mat_t),ncol=length(setdiff(useful_gene_comb_pro2022,colnames(pro_2022_valid_mat_t))))
colnames(add_mat) = setdiff(useful_gene_comb_pro2022,colnames(pro_2022_valid_mat_t))
pro_2022_valid_mat_t_add = (cbind(pro_2022_valid_mat_t,add_mat))




df_cal_2 = merge(pro2022_valid_survival, 
                 data.frame(predict(res_list_comb_pro2022_no_divi[[23]]$m1,newx=pro_2022_valid_mat_t_add[,useful_gene_comb_pro2022],type="response"))%>%
                     rownames_to_column('sample'))


surv_res = plot_survival_contin(df_cal_2, sprintf('s%s',10), plot=T, title=NULL, OS_RF='OS')
surv_res$plot



plot_survival_contin_with_cutpoint(df_cal_2, 
                                   's10', cutpoint='median', cutpoint_type='median', 
                                   plot=T, title=NULL, OS_RF='OS')$plot+
    labs(title=sprintf('%s_%s',23,'s10'))

plot_survival_contin_with_cutpoint(df_cal_2, 
                                   's10', cutpoint=0.3, cutpoint_type='ahead_back', 
                                   plot=T, title=NULL, OS_RF='OS')$plot+
    labs(title=sprintf('%s_%s',23,'s10'))



###################
###################
# valid tcga
###################
###################


exp_mat_tcga = read_tsv('/home/lgh/my_project/data/TCGA/LAML/TCGA-LAML.htseq_fpkm_cleaned.tsv')%>%column_to_rownames('Hugo_Symbol')
exp_mat_tcga = scale(t(exp_mat_tcga))
surv_tcga_df = read_tsv('/home/lgh/my_project/data/TCGA/LAML/TCGA-LAML.survival.cleaned.tsv')
surv_tcga_df[[1]] = str_replace_all(surv_tcga_df[[1]],'[^[:alnum:]]','.')

add_mat = matrix(0,nrow=nrow(exp_mat_tcga),ncol=length(setdiff(useful_gene_comb_pro2022,colnames(exp_mat_tcga))))
colnames(add_mat) = setdiff(useful_gene_comb_pro2022,colnames(exp_mat_tcga))
exp_mat_tcga_add = (cbind(exp_mat_tcga,add_mat))



df_cal_2 = merge(surv_tcga_df, 
                 data.frame(predict(res_list_comb_pro2022_no_divi[[23]]$m1,newx=exp_mat_tcga_add[,useful_gene_comb_pro2022],type="response"))%>%
                     rownames_to_column('sample'))

surv_res = plot_survival_contin(df_cal_2, sprintf('s%s',10), plot=T, title=NULL, OS_RF='OS')
surv_res$plot

plot_survival_contin_with_cutpoint(df_cal_2, 
                                   's10', cutpoint='median', cutpoint_type='median', 
                                   plot=T, title=NULL, OS_RF='OS')$plot+
    labs(title=sprintf('%s_%s',23,'s10'))

plot_survival_contin_with_cutpoint(df_cal_2, 
                                   's10', cutpoint=0.2, cutpoint_type='ahead_back', 
                                   plot=T, title=NULL, OS_RF='OS')$plot














