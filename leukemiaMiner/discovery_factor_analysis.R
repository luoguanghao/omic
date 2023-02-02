library(FactoMineR)
library(factoextra)
library(rvest)
library(tidyverse)


discovery.PCA.factor_analysis <- function(exp_mat_cc,ncp=30, fig=FALSE){

    ###主成分分析

    la_pca<-exp_mat_cc%>%
                PCA(quanti.sup=7:10, ncp = ncp, scale.unit=TRUE)

    var <- get_pca_var(la_pca)

    fviz_pca_var_fig = fviz_pca_var(la_pca,col.var="contrib",
                    gradient.cols=c("#00AFBB","#E7B800","#FC4E07"),
                    repel=TRUE)

    ind <- get_pca_ind(la_pca)

    fviz_pca_ind_cos_fig = fviz_pca_ind(la_pca,col.ind="cos2",select.ind=list(cos2=30),
                        gradient.cols=c("#00AFBB","#E7B800","#FC4E07"), repel=TRUE)

    fviz_pca_ind_fig = fviz_pca_ind(la_pca)

    fviz_screeplot_fig = fviz_screeplot(la_pca,addlabels=TRUE)

    result = list(
        la_pca=la_pca,
        fviz_pca_ind_cos_fig=fviz_pca_ind_cos_fig,
        fviz_pca_ind_fig=fviz_pca_ind_fig,
        fviz_screeplot_fig=fviz_screeplot_fig,
        pca_df = data.frame(get_pca_ind(la_pca)$coord)%>%rownames_to_column('sample')
    )

    discovery_result = list(
        method='PCA',
        result=result,
        factor_df=result$pca_df
    )

    return(discovery_result)


}










