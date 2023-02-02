
# continuous variable
library(glmnet)

#
# disc_y is a data.frame with two column:sample,value
#
pred.elastic_net_cv <- function(disc_mat, disc_y, pred_mat_list, surv_list=NULL, seed=NULL,quiet=FALSE, s='lambda.min',divid_percent=0.7){

    colnames(disc_y) = c('sample','value')

    if(is.null(names(pred_mat_list))){
        names(pred_mat_list) = 1:length(pred_mat_list)
    }

    if(!is.null(seed)){
        set.seed(seed)
        # message(sprintf('set seed=%s',seed))
    }

    colnames(disc_mat) = str_replace_all(colnames(disc_mat),'[^[:alnum:]]','.')
    for(ii in 1:length(pred_mat_list)){
        colnames(pred_mat_list[[ii]]) = str_replace_all(colnames(pred_mat_list[[ii]]),'[^[:alnum:]]','.')
    }    

    # 使得样本顺序一样
    disc_mat = disc_mat[intersect(rownames(disc_mat),disc_y$sample),]
    disc_y = (disc_y%>%column_to_rownames('sample'))[intersect(rownames(disc_mat),disc_y$sample),]
    # ===

    useful_genes = colnames(disc_mat)
    for(ii in 1:length(pred_mat_list)){
        useful_genes = intersect( useful_genes,colnames(pred_mat_list[[ii]]) )
    }

    sample_res = sample(1:length(disc_y),length(disc_y)*(1-divid_percent),replace = FALSE)
    #print('D')
    if(divid_percent<1){
        x_valid = disc_mat[sample_res,useful_genes]
        y_valid = disc_y[sample_res]
    }
    #print('C')
    x_do = disc_mat[setdiff(1:length(disc_y),sample_res),useful_genes]
    y_do = disc_y[setdiff(1:length(disc_y),sample_res)]

    cvob1 = cv.glmnet(x_do, y_do)
    if(!quiet) plot(cvob1)

    do_res = (cor.test(
        (stats::predict(cvob1, newx = x_do, s = s))[,1],y_do
    ))
    #print('B')
    valid_res = NULL
    if(divid_percent<1){
        valid_res = (cor.test(
            (stats::predict(cvob1, newx = x_valid, s = s))[,1],y_valid
        ))
    }
    #print('A')
    if(!quiet){
        print(do_res)
        plot(
            (stats::predict(cvob1, newx = x_do, s = s))[,1],y_do
        )

        if(divid_percent<1){
            print(valid_res)
            plot(
                (stats::predict(cvob1, newx = x_valid, s = s))[,1],y_valid
            )
        }
    }
    #print('E')
    pred.result.list = list()
    for(ii in 1:length(pred_mat_list)){

        tmp_predict_res = stats::predict(cvob1, newx = as.matrix(pred_mat_list[[ii]])[,useful_genes], s = s)
        pred.result.list[[names(pred_mat_list)[[ii]]]] = data.frame(tmp_predict_res)%>%rownames_to_column('sample')
        colnames(pred.result.list[[names(pred_mat_list)[[ii]]]])[2]  = 'pred_value'

    }

    fit <- glmnet(x_do, y_do)

    return(list(cvob1=cvob1, 
                train_cor_test=list(do_res,do_res,valid_res=valid_res), 
                result=pred.result.list, 
                useful_genes=useful_genes,fit=fit)
        )


}

# batch_survival_after_factor_analysis



if(FALSE){

    for(seed_tmp in 501:1000){
        pred.elastic_net_cv_res = pred.elastic_net_cv(disc_mat=t(pro_mat[(sel_pro_by_na),factor_analysis_res_N$factor_df$sample]), 
                            disc_y=factor_analysis_res_N$factor_df$Dim.10, 
                            pred_mat_list=list(AML2022=pro_2022_mat_t,AML2022_valid=pro_2022_valid_mat_t), seed=seed_tmp, quiet=TRUE)
        
        df_for_cal1 = merge(pred.elastic_net_cv_res$result$AML2022,pro2022_survival,by='sample')
        df_for_cal2 = merge(pred.elastic_net_cv_res$result$AML2022_valid,pro2022_valid_survival,by='sample')
        if(plot_survival_contin_with_cutpoint(df_for_cal1, sprintf('pred_value'),cutpoint='median', cutpoint_type='median', plot=F, title=NULL, OS_RF='OS')$pval<0.05 &
        plot_survival_contin_with_cutpoint(df_for_cal2, sprintf('pred_value'),cutpoint='median', cutpoint_type='median', plot=T, title=NULL, OS_RF='OS')$pval<0.05
        ){
            # print(seed_tmp)
            df_for_cal = merge(pred.elastic_net_cv_res$result$AML2022,pro2022_survival,by='sample')
            print(plot_survival_contin_with_cutpoint(df_for_cal, sprintf('pred_value'),cutpoint='median', cutpoint_type='median', plot=T, title=NULL, OS_RF='OS')$plot)
            df_for_cal = merge(pred.elastic_net_cv_res$result$AML2022_valid,pro2022_valid_survival,by='sample')
            print(plot_survival_contin_with_cutpoint(df_for_cal, sprintf('pred_value'),cutpoint='median', cutpoint_type='median', plot=T, title=NULL, OS_RF='OS')$plot)
            
        }
        
    }

}


## elastic net coxph
#
# df_cal_1 : sample,OS,OSS,s0,s1~sn
# df_cal_2_list
# res_df_1 : logrank p value: 'ldb''p''n_high''n_low'
# res_df_2_list
# res_df_coxph_1 : coxph result: 'beta''HR..95..CI.for.HR.''wald.test''p.value'
# res_df_coxph_2_list
# surv_res_1 : logRank生存分析的原始结果
# surv_res_2_list
#
library(fastcox)
pred.elastic_net_coxph_cv <- function(disc_mat, disc_surv, pred_mat_list, pred_surv_list, seed=NULL,quiet=FALSE,OS_RF='OS'){

    if(is.null(names(pred_mat_list))){
        names(pred_mat_list) = 1:length(pred_mat_list)
    }

    if(!is.null(seed)){
        set.seed(seed)
        # message(sprintf('set seed=%s',seed))
    }

    useful_genes = colnames(disc_mat)
    for(ii in 1:length(pred_mat_list)){
        useful_genes = intersect( useful_genes,colnames(pred_mat_list[[ii]]) )
    }

    x = disc_mat
    disc_surv = disc_surv%>%filter(sample%in%rownames(disc_mat))
    x = x[disc_surv$sample,]
    if(OS_RF=='OS'){
        y = disc_surv[['OS']]
        status = disc_surv[['OSS']]
    }else{
        y = disc_surv[['RFS']]
        status = disc_surv[['RFSS']]        
    }

    # useful_genes = intersect(colnames(pro_2022_mat_t),colnames(x))
    sample_res = rownames(x)[1:length(y)]

    x_do = x[order(rownames(x)), useful_genes]
    y_do = y[order(rownames(x))]
    status_do = status[order(rownames(x))]

    cv1<-cv.cocktail(x=x_do,y=y_do,d=status_do,alpha=0.5,nfolds=10) 

    lbd_seq = seq(cv1$lambda.min, cv1$lambda.1se, (cv1$lambda.1se-cv1$lambda.min)/20)

    print(c('##',lbd_seq))
    #m1 = list()
    if(length(lbd_seq)!=1){
    
        m1<-cocktail(x=x_do,y=y_do,d=status_do,alpha=0.5,lambda=lbd_seq)

        #selected_featrue = as.matrix(predict(m1,newx=x_do,type="coefficients"))[,'s40']
        #selected_featrue = selected_featrue[selected_featrue>0]

        res_df_1 = list(ldb=c(),p=c(),n_high=c(),n_low=c())
        res_df_coxph_1 = list(beta=c(),`HR (95% CI for HR)`=c(),wald.test=c(),p.value=c())
        df_cal_1 = merge(disc_surv,data.frame(predict(m1,newx=x_do,type="response"))%>%rownames_to_column('sample'))
        surv_res_1 = list()
        for(i in 1:(length(lbd_seq)-1)) {
            surv_res = plot_survival_contin(df_cal_1, sprintf('s%s',i), plot=F, title=NULL, OS_RF=OS_RF)

            res_df_1$ldb = c(res_df_1$ldb, sprintf('s%s',i))
            res_df_1$p = c(res_df_1$p, surv_res$pval)

            res_df_1$n_high = c(res_df_1$n_high, surv_res$surv_diff$n[[1]])
            res_df_1$n_low = c(res_df_1$n_low, surv_res$surv_diff$n[[2]])

            surv_res_1[[i]] = surv_res
            #print(plot_survival_contin(df_cal, sprintf('s%s',i), plot=T, title=NULL, OS_RF='OS')$plot+labs(title=sprintf('s%s',i)))


            res_df_coxph_1_tmp =  summary(coxph(as.formula(sprintf('Surv(%s,%sS)~s%s',OS_RF,OS_RF,i)), df_cal_1))
            p.value<-signif(res_df_coxph_1_tmp$wald["pvalue"], digits=2)
            wald.test<-signif(res_df_coxph_1_tmp$wald["test"], digits=2)
            beta<-signif(res_df_coxph_1_tmp$coef[1], digits=2);#coeficient beta
            HR <-signif(res_df_coxph_1_tmp$coef[2], digits=2);#exp(beta)
            HR.confint.lower <- signif(res_df_coxph_1_tmp$conf.int[,"lower .95"],2)
            HR.confint.upper <- signif(res_df_coxph_1_tmp$conf.int[,"upper .95"],2)
            HR <- paste0(HR, " (", 
                        HR.confint.lower, "-", HR.confint.upper, ")")

            res_df_coxph_1[['beta']] = c(res_df_coxph_1[['beta']],beta)
            res_df_coxph_1[['HR (95% CI for HR)']] = c(res_df_coxph_1[['HR (95% CI for HR)']],HR)
            res_df_coxph_1[['wald.test']] = c(res_df_coxph_1[['wald.test']],wald.test)
            res_df_coxph_1[['p.value']] = c(res_df_coxph_1[['p.value']],p.value)



        }
        res_df_1 = data.frame(res_df_1)
        res_df_coxph_1 = data.frame(res_df_coxph_1)

        # ==========================
        # ==========================

        df_cal_2_list = list()
        res_df_2_list = list()
        res_df_coxph_2_list = list()
        surv_res_2_list = list()
        for(ii in 1:length(pred_mat_list)){

            res_df_2 = list(ldb=c(),p=c(),n_high=c(),n_low=c())
            res_df_coxph_2 = list(beta=c(),`HR (95% CI for HR)`=c(),wald.test=c(),p.value=c())
            df_cal_2 = merge(pred_surv_list[[ii]], 
                            data.frame(predict(m1,newx=pred_mat_list[[ii]][,useful_genes],type="response"))%>%
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
                    surv_res = plot_survival_contin(df_cal_2, sprintf('s%s',i), plot=F, title=NULL, OS_RF=OS_RF)
                    #print(plot_survival_contin(df_cal, sprintf('s%s',i), plot=T, title=NULL, OS_RF='OS')$plot+labs(title=sprintf('s%s',i)))
                    res_df_2$ldb = c(res_df_2$ldb, sprintf('s%s',i))
                    res_df_2$p = c(res_df_2$p, surv_res$pval)

                    res_df_2$n_high = c(res_df_2$n_high, surv_res$surv_diff$n[[1]])
                    res_df_2$n_low = c(res_df_2$n_low, surv_res$surv_diff$n[[2]])

                    surv_res_2[[i]] = surv_res
                }

                res_df_coxph_2_tmp =  summary(coxph(as.formula(sprintf('Surv(%s,%sS)~s%s',OS_RF,OS_RF,i)), df_cal_2))
                p.value<-signif(res_df_coxph_2_tmp$wald["pvalue"], digits=2)
                wald.test<-signif(res_df_coxph_2_tmp$wald["test"], digits=2)
                beta<-signif(res_df_coxph_2_tmp$coef[1], digits=2);#coeficient beta
                HR <-signif(res_df_coxph_2_tmp$coef[2], digits=2);#exp(beta)
                HR.confint.lower <- signif(res_df_coxph_2_tmp$conf.int[,"lower .95"],2)
                HR.confint.upper <- signif(res_df_coxph_2_tmp$conf.int[,"upper .95"],2)
                HR <- paste0(HR, " (", 
                            HR.confint.lower, "-", HR.confint.upper, ")")

                res_df_coxph_2[['beta']] = c(res_df_coxph_2[['beta']],beta)
                res_df_coxph_2[['HR (95% CI for HR)']] = c(res_df_coxph_2[['HR (95% CI for HR)']],HR)
                res_df_coxph_2[['wald.test']] = c(res_df_coxph_2[['wald.test']],wald.test)
                res_df_coxph_2[['p.value']] = c(res_df_coxph_2[['p.value']],p.value)

            }
            res_df_2 = data.frame(res_df_2)
            res_df_coxph_2 = data.frame(res_df_coxph_2)

            #res_list_comb_pro2022[[ii]]$df_cal_pro2022 = df_cal_2
            #res_list_comb_pro2022[[ii]]$res_df_pro2022 = res_df_2
            #res_list_comb_pro2022[[ii]]$surv_res_pro2022 = surv_res_2

            df_cal_2_list[[ names(pred_mat_list)[[ii]] ]] = df_cal_2
            res_df_2_list[[ names(pred_mat_list)[[ii]] ]] = res_df_2
            res_df_coxph_2_list[[ names(pred_mat_list)[[ii]] ]] = res_df_coxph_2
            surv_res_2_list[[ names(pred_mat_list)[[ii]] ]] = surv_res_2

        }

        tmp_res = list(
            sample_res=sample_res,
            cv1=cv1,
            m1=m1,
            #selected_featrue=selected_featrue,

            x_do=x_do,
            y_do=y_do,
            status_do=status_do,

            df_cal_1=df_cal_1,
            df_cal_2_list=df_cal_2_list,

            res_df_1 = res_df_1,
            res_df_2_list = res_df_2_list,

            res_df_coxph_1 = res_df_coxph_1,
            res_df_coxph_2_list = res_df_coxph_2_list,

            surv_res_1 = surv_res_1,
            surv_res_2_list = surv_res_2_list
        )
        print('%%%%%%%$$')
    }else{
        tmp_res = list(
            sample_res=sample_res,
            cv1=cv1,

            x_do=x_do,
            y_do=y_do,
            status_do=status_do,
            # m1=m1,
            #selected_featrue=selected_featrue,

            df_cal_1="NO",
            df_cal_2_list="NO",

            res_df_1 = "NO",
            res_df_2_list = "NO",

            res_df_coxph_1 = "NO",
            res_df_coxph_2_list = "NO",

            surv_res_1 = "NO",
            surv_res_2_list = "NO"
        )        
    }

    return(tmp_res)



}








# pred.elastic_net_cv

pred.randomForest_cv <- function(disc_mat, disc_group_df, pred_mat_list, pred_surv_list=NULL, use_feature_n=NULL, seed=NULL,quiet=FALSE,OS_RF='OS',divid_percent=0.8){

    # group_df = result_list_HCPC[[102]][[3]]$res_post_subtyping$group_df
    # disc_mat = scale(pro_mat_t[,sel_pro_by_na])
    # pred_mat_list = list(pro_2022_mat_t=pro_2022_mat_t, 
    #                     pro_2022_valid_mat_t=pro_2022_valid_mat_t)
    # seed = 23

    colnames(disc_mat) = str_replace_all(colnames(disc_mat),'[^[:alnum:]]','.')
    for(ii in 1:length(pred_mat_list)){
        colnames(pred_mat_list[[ii]]) = str_replace_all(colnames(pred_mat_list[[ii]]),'[^[:alnum:]]','.')
    }  

    useful_genes = colnames(disc_mat)
    for(ii in 1:length(pred_mat_list)){
        useful_genes = intersect( useful_genes,colnames(pred_mat_list[[ii]]) )
    }

    mat_for_ml = merge(data.frame(disc_mat[,useful_genes])%>%rownames_to_column('sample'),disc_group_df )


    #将总数据集分为训练集（占 70%）和测试集（占 30%）

    set.seed(seed) # <-------------

    mat_for_ml = mat_for_ml%>%column_to_rownames('sample')
    mat_for_ml$group = as.factor(mat_for_ml$group)

    if(divid_percent==1){
        
        mat_for_train <- mat_for_ml

        train.forest <- randomForest::randomForest(group ~ ., data = mat_for_train, importance = TRUE, ntree=500)

        # train.forest
        print(train.forest)
        message('======================')
        train_predict <- predict(train.forest, mat_for_train)
        compare_train <- table(train_predict, mat_for_train$group)
        print(compare_train)
        # sum(diag(compare_train)/sum(compare_train))    
        message('======================')
    }else{
        select_train <- sample(nrow(mat_for_ml), nrow(mat_for_ml)*0.8)
        mat_for_train <- mat_for_ml[select_train, -1]
        mat_for_test <- mat_for_ml[-select_train, -1]

        train.forest <- randomForest::randomForest(group ~ ., data = mat_for_train, importance = TRUE, ntree=500)

        # train.forest
        print(train.forest)
        message('======================')
        train_predict <- predict(train.forest, mat_for_train)
        compare_train <- table(train_predict, mat_for_train$group)
        print(compare_train)
        # sum(diag(compare_train)/sum(compare_train))    
        message('======================')
        #使用测试集评估
        test_predict <- predict(train.forest, mat_for_test)
        compare_test <- table(mat_for_test$group, test_predict, dnn = c('Actual', 'Predicted'))
        print(compare_test)
        message('======================\n======================')

        #result_list[[i_seed]][['train.forest']] = train.forest
        #result_list[[i_seed]][['train_predict']] = train_predict
        #result_list[[i_seed]][['compare_train']] = compare_train
        #result_list[[i_seed]][['test_predict']] = test_predict
        #result_list[[i_seed]][['compare_test']] = compare_test
    }

    ####################
    #################### predict
    # pro2022_survival pro2022_valid_survival
    predict_res_list = list()

    for(i in 1:length(pred_mat_list)){

        predict_res <- predict(train.forest, pred_mat_list[[i]][,useful_genes])
        predict_res = data.frame(predict_res)%>%rownames_to_column('sample')
        colnames(predict_res)[2] = 'group'

        predict_res$group = as.character(predict_res$group)

        predict_res_list[[names(pred_mat_list)[i]]] = predict_res

    }

    return(list(predict_res_list=predict_res_list,train.forest=train.forest,useful_genes=useful_genes))

}

























