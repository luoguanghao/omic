library("survival")
library("survminer")
library(ggplot2)



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


plot_survival_binary <- function(df_for_suv_plot, var, plot=FALSE){

    f <- as.formula(paste(c('Surv(OS, OSS)',var),collapse=' ~ '))
    fit <- surv_fit(f, data = df_for_suv_plot)
    #fit <- surv_fit(Surv(df_for_suv_plot[['OS']], df_for_suv_plot[['OSS']]) ~ df_for_suv_plot[[var]])

    if(plot){       
        ggsurvplot(fit,
            pval = TRUE, conf.int = TRUE,
            risk.table = TRUE, # Add risk table
            risk.table.col = "strata", # Change risk table color by groups
            linetype = "strata", # Change line type by groups
            surv.median.line = "hv", # Specify median survival
            ggtheme = theme_bw(), # Change ggplot2 theme
            palette = c("#E7B800", "#2E9FDF")
            )
    }

    #surv_diff <- survdiff(Surv(df_for_suv_plot[['OS']], df_for_suv_plot[['OSS']]) ~ df_for_suv_plot[[var]])
    #pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

    #return(list(surv_diff=surv_diff, pval=pval, surv_fit=fit))
}


# 生存分析 连续变量，寻找最佳cutoff #
## df_for_suv_plot: df: sample;some binary var; OS; OSS
## var: str, descript the col name of continue var
#
# df_for_suv_plot = na.omit(merge(ssgsea_senunsen_df, suv_df))
plot_survival_contin <- function(df_for_suv_plot, var, plot=FALSE){

    res.cut <- surv_cutpoint(df_for_suv_plot, time = "OS", event = "OSS",
    variables = c(var)) ; res.cut

    res.cat <- surv_categorize(res.cut)

    f <- as.formula(paste(c('Surv(OS, OSS)',var),collapse=' ~ '))

    fit <- surv_fit(f, data = res.cat)

    if(plot){
        ggsurvplot(fit,
            pval = TRUE, conf.int = TRUE,
            risk.table = TRUE, # Add risk table
            risk.table.col = "strata", # Change risk table color by groups
            linetype = "strata", # Change line type by groups
            surv.median.line = "hv", # Specify median survival
            ggtheme = theme_bw(), # Change ggplot2 theme
            palette = c("#E7B800", "#2E9FDF"))
    }

    surv_diff <- survdiff(f, data = lung)
    pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

    return(list(surv_diff=surv_diff, pval=pval, surv_fit=fit))

}



