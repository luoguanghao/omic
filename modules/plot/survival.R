library("survival")
library("survminer")
library(ggplot2)

library(gridExtra)
library("gtable")
library("grid")
library(ggthemes)

######### 这两个函数不能用，问题出在get(var)那里


# only survival analysis without plot
if(FALSE){
    surv_diff <- survdiff(Surv(df_for_suv_plot[['OS']], df_for_suv_plot[['OSS']]) ~ df_for_suv_plot[[var]])
    pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)    
}

# 生存分析 二元变量 #
## df_for_suv_plot: df: sample;some binary var; OS; OSS
## var: str, descript the col name of binary var
#
# df_for_suv_plot = na.omit(merge(ssgsea_senunsen_df, suv_df))


plot_survival_binary <- function(df_for_suv_plot, var, plot=FALSE, conf.int=FALSE, OS_RF='OS'){
    if(OS_RF=='OS'){
        f <- as.formula(paste(c('Surv(OS, OSS)',var),collapse=' ~ '))
    }else{
        f <- as.formula(paste(c('Surv(RFS, RFSS)',var),collapse=' ~ '))
    }
    fit <- surv_fit(f, data = df_for_suv_plot)
    #fit <- surv_fit(Surv(df_for_suv_plot[['OS']], df_for_suv_plot[['OSS']]) ~ df_for_suv_plot[[var]])

    surv_diff <- survdiff(f, data = df_for_suv_plot)
    pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

    # ==
    surv_diff_pairwise <- pairwise_survdiff(f, data = df_for_suv_plot)

    df_plot = data.frame(sapply(as.data.frame(surv_diff_pairwise$p.value), function(x) signif(x,3)))
    rownames(df_plot) = rownames(surv_diff_pairwise$p.value)
    df_plot = df_plot%>%rownames_to_column('.')

    # == cox

    cox_res = coxph(f, data=df_for_suv_plot)

    beta<-signif(summary(cox_res)$coef[1], digits=2)
    HR <-signif(summary(cox_res)$coef[2], digits=2)
    lower = signif(summary(cox_res)$conf.int[,"lower .95"],2)
    upper = signif(summary(cox_res)$conf.int[,"upper .95"],2)
    p.value<-signif(summary(cox_res)$wald["pvalue"], digits=2)


    p=NULL
    pp=NULL
    if(plot){       
        p=ggsurvplot(fit,
            pval = TRUE, conf.int = conf.int,
            risk.table = TRUE, # Add risk table
            risk.table.col = "strata", # Change risk table color by groups
            linetype = "strata", # Change line type by groups
            # surv.median.line = "hv", # Specify median survival
            #ggtheme = theme_bw(), # Change ggplot2 theme
            palette = c( "#FFCC00","#FF0033")
            #palette = c("cell")
            )


        df_plot.p <- ggtexttable(df_plot, rows = NULL, 
                        theme = ttheme("mOrange"))

        pp <- arrangeGrob(p[[1]], p[[2]]+theme(legend.position="none"), df_plot.p+theme_void(), 
                    ncol = 1, nrow = 5, 
                    layout_matrix = rbind(c(1),c(1), c(2),c(3)))


    }


    return(list(surv_diff=surv_diff, pval=pval, surv_fit=fit, df_pairwise=df_plot, plot=p, plot_with_df=pp,
                df_for_suv_plot=df_for_suv_plot, cox_res = list(beta=beta,HR=HR,lower=lower,upper=upper,p.value=p.value)))
}


# 生存分析 连续变量，寻找最佳cutoff #
## df_for_suv_plot: df: sample;some binary var; OS; OSS
## var: str, descript the col name of continue var
#
# df_for_suv_plot = na.omit(merge(ssgsea_senunsen_df, suv_df))
plot_survival_contin <- function(df_for_suv_plot, var, plot=FALSE, title=NULL, OS_RF='OS'){

    # res.cut <- surv_cutpoint(df_for_suv_plot, time = "OS", event = "OSS",

    if(OS_RF=='OS'){
        res.cut <- surv_cutpoint(df_for_suv_plot, time = "OS", event = "OSS", variables = c(var)) ; res.cut
    }else{
        res.cut <- surv_cutpoint(df_for_suv_plot, time = "RFS", event = "RFSS", variables = c(var)) ; res.cut
    }

    

    res.cat <- surv_categorize(res.cut)

    # f <- as.formula(paste(c('Surv(OS, OSS)',var),collapse=' ~ '))

    if(OS_RF=='OS'){
        f <- as.formula(paste(c('Surv(OS, OSS)',var),collapse=' ~ '))
    }else{
        f <- as.formula(paste(c('Surv(RFS, RFSS)',var),collapse=' ~ '))
    }

    fit <- surv_fit(f, data = res.cat)
    p=NULL
    if(plot){
        p=ggsurvplot(fit,
            pval = TRUE, conf.int = FALSE,
            risk.table = TRUE, # Add risk table
            risk.table.col = "strata", # Change risk table color by groups
            linetype = "strata", # Change line type by groups
            # surv.median.line = "hv", # Specify median survival
            ggtheme = theme_few(), # Change ggplot2 theme
            #palette = c( "#FFCC00","#003399")
            palette = c("cell") 
            )
        # print(p)
    }

    surv_diff <- survdiff(f, data = res.cat)
    pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

    return(list(surv_diff=surv_diff, pval=pval, surv_fit=fit,res.cut=res.cut,res.cat=res.cat,plot=p))

}

######### 
# < cutpoint is low
# ahead_back smaller than 0.5
# cutpoint = value
plot_survival_contin_with_cutpoint <- function(df_for_suv_plot, var, cutpoint='median', cutpoint_type='median', plot=FALSE, title=NULL, OS_RF='OS'){

    if(OS_RF=='OS'){

        if(cutpoint_type=='ahead_back'){

            cutpoint1 = quantile(df_for_suv_plot[[var]], cutpoint)
            cutpoint2 = quantile(df_for_suv_plot[[var]], 1-cutpoint)

            tmp = df_for_suv_plot[,c('sample','OSS','OS',var)]
            tmp$group = 'mid'

            tmp[tmp[[var]] <= cutpoint1,]$group = 'low'      
            tmp[tmp[[var]] >= cutpoint2,]$group = 'high'   
            tmp = tmp%>%filter(group!='mid')
            res = plot_survival_binary(tmp, 'group', plot=T, conf.int=FALSE, OS_RF='OS')

            return(res)


        }

        if(cutpoint_type=='median'){
            cutpoint = median(df_for_suv_plot[[var]])
        }else if(cutpoint_type=='quantile'){
            cutpoint = quantile(df_for_suv_plot[[var]], cutpoint)
        }

        tmp = df_for_suv_plot[,c('sample','OSS','OS',var)]
        tmp$group = 'high'

        tmp[tmp[[var]] < cutpoint,]$group = 'low'

        res = plot_survival_binary(tmp, 'group', plot=T, conf.int=FALSE, OS_RF='OS')

        return(res)
    }

    if(OS_RF=='RFS'){

        if(cutpoint_type=='ahead_back'){

            cutpoint1 = quantile(df_for_suv_plot[[var]], cutpoint)
            cutpoint2 = quantile(df_for_suv_plot[[var]], 1-cutpoint)

            tmp = df_for_suv_plot[,c('sample','RFSS','RFS',var)]
            tmp$group = 'mid'

            tmp[tmp[[var]] <= cutpoint1,]$group = 'low'      
            tmp[tmp[[var]] >= cutpoint2,]$group = 'high'   
            tmp = tmp%>%filter(group!='mid')
            res = plot_survival_binary(tmp, 'group', plot=T, conf.int=FALSE, OS_RF='RFS')

            return(res)


        }

        if(cutpoint_type=='median'){
            cutpoint = median(df_for_suv_plot[[var]])
        }else if(cutpoint_type=='quantile'){
            cutpoint = quantile(df_for_suv_plot[[var]], cutpoint)
        }

        tmp = df_for_suv_plot[,c('sample','RFSS','RFS',var)]
        tmp$group = 'high'

        tmp[tmp[[var]] < cutpoint,]$group = 'low'

        res = plot_survival_binary(tmp, 'group', plot=T, conf.int=FALSE, OS_RF='RFS')

        return(res)
    }

}







