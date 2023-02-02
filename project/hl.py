





def do_elastic_net(X, y, gene_ls, alphas, cv, max_iter, random_state, selection, verbose=1):
    '''
    X: sample*gene np.array
    y: sample np.array
    gene_ls:np.array
    '''
    regr = ElasticNetCV(alphas=alphas, cv=cv, max_iter=max_iter, random_state=random_state, selection=selection)
    # print(X.shape)
    regr.fit(X, y)

    # predict_y = regr.predict(X)
    '''
    feature_ls = []
    feature_weight_df = []
    for i in range(len(regr.coef_)):
        if regr.coef_[i]!=0:
            feature_ls.append(all_feature_ls[i])
            feature_weight_df.append([all_feature_ls[i], regr.coef_[i]])
    feature_weight_df = pd.DataFrame(feature_weight_df,columns=['Feature','Weight']).set_index(['Feature'])
    '''
    feature_weight_df = pd.DataFrame({'Feature':gene_ls[regr.coef_!=0],
                                    'Weight':regr.coef_[regr.coef_!=0]}).set_index(['Feature'])
    feature_ls = np.array(feature_weight_df.index)


    if verbose:  
        print('Elastic Net keep %s features...'%len(feature_ls))
    
    return {'model':regr, 'feature_ls':feature_ls, 
            'feature_weight_df':feature_weight_df, 
            'X':X, 'y':y}


def plot_post_feature_selection(expr_mat, drug_response, importance_feature):
    mat_for_plot = expr_mat.T.loc[:,importance_feature].iloc[140:190]

    %R library(pheatmap)
    %R -i mat_for_plot
    %R p_res = pheatmap(mat_for_plot)

    order = %R p_res$tree_row$order
    labels = %R p_res$tree_row$label
    plt.plot(all_drug_sen_df.loc[drug].set_index('lab_id').loc[drug_expr_sample_ls,'auc'].loc[labels[order-1]])


def do_simple_en(expr_mat, drug_response, log_flag, 
                en_max_iter, feature_number_threshold, verbose=1):
    '''
    expr_mat   row: gene  columns:sample
    drug_response  Series
    '''
    drug_sample_ls = list(drug_response.index)
    expr_sample_ls = list(expr_mat.columns)
    drug_expr_sample_ls = list( set(drug_sample_ls)&set(expr_sample_ls) )
    
    expr_mat = expr_mat[drug_expr_sample_ls]
    drug_response = drug_response.loc[drug_expr_sample_ls]

    gene_list = np.array(expr_mat.index)


    if log_flag is 1:
        X_train = np.log( np.array(expr_mat.T)+1 )
    else:
        X_train = np.array(expr_mat.T)

    y_train = np.array(drug_response)

    
    #########
    # model #
    #########
    en_res = do_elastic_net(X=X_train, y=y_train, gene_ls=gene_list, alphas=None, cv=10, 
                            random_state=random.randint(1,100), selection='random', max_iter=en_max_iter )    

    return {'en_res':en_res,
            'drug_expr_sample_ls':drug_expr_sample_ls,
            'gene_list':gene_list}
    ##########
    # plot_post_feature_selection
    ##########
    importance_feature = abs(en_res['feature_weight_df']).sort_values('Weight').index[-feature_number_threshold:]


    mat_for_plot = expr_mat.T.loc[:,importance_feature].iloc[140:190]

    %R library(pheatmap)
    %R -i mat_for_plot
    %R p_res = pheatmap(mat_for_plot)

    order = %R p_res$tree_row$order
    labels = %R p_res$tree_row$label
    plt.plot(all_drug_sen_df.loc[drug].set_index('lab_id').loc[drug_expr_sample_ls,'auc'].loc[labels[order-1]])













