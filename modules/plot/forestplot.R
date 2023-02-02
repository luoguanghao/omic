
if(FALSE){
    df_for_plot = pair_compare_maf_res$pt.vs.rt$result%>%filter(or!=Inf)


    cochrane_from_rmeta <- 
    structure(list(
        mean  = c(NA, NA, df_for_plot$or[1:10]), 
        lower = c(NA, NA, df_for_plot$ci.low[1:10]),
        upper = c(NA, NA, df_for_plot$ci.up[1:10])),
        .Names = c("mean", "lower", "upper"), 
        row.names = c(NA, -12L), 
        class = "data.frame")

    tabletext<-cbind(
    c("", "Study", df_for_plot[['Hugo_Symbol']][1:10]),
    c("< -1", "(< -1)", df_for_plot[['< -1']][1:10]),
    c("> -0.7	", "(> -0.7	)", df_for_plot[['> -0.7']][1:10]),
    c("", "OR",sprintf(

    forestplot(tabletext, 
            cochrane_from_rmeta,new_page = TRUE,
            is.summary=c(TRUE,TRUE,rep(FALSE,8),TRUE),
            clip=c(0.1,2.5), 
            xlog=TRUE, 
            col=fpColors(box="royalblue",line="darkblue", summary="royalblue"))





}



