


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
                                    'Weight':regr.coef_[regr.coef_!=0]}).set_index(['Feature']).sort_values('Weight')
    feature_ls = np.array(feature_weight_df.index)


    if verbose:  
        print('Elastic Net keep %s features...'%len(feature_ls))
    
    return {'model':regr, 'feature_ls':feature_ls, 
            'feature_weight_df':feature_weight_df, 
            'X':X, 'y':y}






























