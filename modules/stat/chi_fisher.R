



do_fisher_for_pair <- function(group_df,feature_df,feature){
    feature_df$feature = feature_df[[feature]]
    df_for_cal = merge(group_df,feature_df)
    fisher_matrix = matrix(c(0,0,0,0),ncol=2,
                        dimnames = list(c("this_feature", "no_this_feature"),
                                        c("this_group", "no_this_group"))
                )
    
    fisher_res_df = list( group=c(),feature=c(),p=c(),

                     this_feature_in_this_group=c(),this_feature_in_no_this_group=c() ) 
    for(g1 in unique(df_for_cal$group)){

        for(f1 in unique(df_for_cal$feature)){

            fisher_matrix[1,1] = dim(df_for_cal%>%filter(group==g1&feature==f1))[1]

            fisher_matrix[1,2] = dim(df_for_cal%>%filter(group!=g1&feature==f1))[1]

            fisher_matrix[2,1] = dim(df_for_cal%>%filter(group==g1&feature!=f1))[1]

            fisher_matrix[2,2] = dim(df_for_cal%>%filter(group!=g1&feature!=f1))[1]

            #fisher_res_df$param = c(fisher_res_df$param,param_ls[param_i])

            #fisher_res_df$grouping_k = c(fisher_res_df$grouping_k,grouping_ls[grouping_i])

            fisher_res_df$group = c(fisher_res_df$group,g1)

            fisher_res_df$feature = c(fisher_res_df$feature,f1)

            fisher_res_df$p = c(fisher_res_df$p,fisher.test(fisher_matrix)$p.value)

            fisher_res_df$this_feature_in_this_group = c(fisher_res_df$this_feature_in_this_group,

                                                            sprintf('%s/%s',fisher_matrix[1,1],fisher_matrix[1,1]+fisher_matrix[2,1])

                                                        )

            fisher_res_df$this_feature_in_no_this_group = c(fisher_res_df$this_feature_in_no_this_group,

                                                            sprintf('%s/%s',fisher_matrix[1,2],fisher_matrix[1,2]+fisher_matrix[2,2])

                                                            )

        }    

    }

    fisher_res_df = data.frame(fisher_res_df)

    fisher_res_df_pairwise = fisher_res_df[order(fisher_res_df$p),]

    return(fisher_res_df_pairwise)

}

do_fisher_overall <- function(group_df,feature_df,feature){
    feature_df$feature = feature_df[[feature]]
    df_for_cal = merge(group_df,feature_df)

    fisher_matrix = matrix(0,

                            length(unique(df_for_cal$feature)),

                            length(unique(df_for_cal$group)),

                            dimnames = list(unique(df_for_cal$feature),

                                        unique(df_for_cal$group))

                            )

        # ===

    for( i_r in 1:length(unique(df_for_cal$feature)) ){

        for( i_c in 1:length(unique(df_for_cal$group)) ){

            fisher_matrix[i_r, i_c] = dim( df_for_cal%>%filter(group==unique(df_for_cal$group)[i_c]&

                                                                feature==unique(df_for_cal$feature)[i_r]) )[1]

        }            

    }

    return( try( fisher.test(fisher_matrix) ) )


}

















