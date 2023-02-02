'''
一些变量的名字要改 比如res_en_2之类的
'''


import itertools
import warnings

import math

import pandas as pd
import pylab
import numpy as np


from sklearn.linear_model import enet_path
from sklearn import preprocessing
from sklearn import model_selection
from sklearn import linear_model # must use the module rather than classes to

import matplotlib.pyplot as plt

import random
%load_ext rpy2.ipython

import pandas as pd
import numpy as np
import math, csv, pickle
import random, os
from sklearn.linear_model import ElasticNetCV
from sklearn.datasets import make_regression

from sklearn.model_selection import train_test_split

import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

###################
## core function ##
###################
def do_elastic_net(X, y, gene_ls, alphas, cv, max_iter, random_state, selection, verbose=1):
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

def train_test_set(df_samlple=None, df_lable=None, sample_ls=None, test_size=0.3, return_type='index'):
    print('return tpye is "%s"'%return_type)
    if return_type=='data':
        sample_ls = df_samlple.index
    
    n_sample = len(sample_ls)
    n_test = round(n_sample*test_size)
    train_index_ls = []
    test_index_ls = []
    train_ls = []
    test_ls = []
    
    test_index_ls = np.array(random.sample(range(n_sample), k=n_test))
    train_index_ls = np.array(list(set(range(n_sample))-set(test_index_ls)))
    train_ls = sample_ls[train_index_ls]
    test_ls = sample_ls[test_index_ls]
    
    if return_type=='index':
        return {'train_index_ls':train_index_ls,
                'test_index_ls':test_index_ls,
                'train_ls':train_ls,
                'test_ls':test_ls}
    
    elif return_type=='data':
        df_samlple_train = df_samlple.loc[train_ls]
        df_samlple_test = df_samlple.loc[test_ls]
        df_label_train = df_lable.loc[train_ls]
        df_label_test = df_lable.loc[test_ls]
        return {'train_index_ls':train_index_ls,
                'test_index_ls':test_index_ls,
                'train_ls':train_ls,
                'test_ls':test_ls,
                'df_samlple_train':df_samlple_train,
                'df_samlple_test':df_samlple_test,
                'df_label_train':df_label_train,
                'df_label_test':df_label_test
               }

def feature_pheatmap(df, feature_ls, tran, df_lable=None):
    '''

    '''
    if tran == 'log':
        genes_samples_mat = np.log(df.loc[:,feature_ls]+0.1)
    else:
        genes_samples_mat = df.loc[:,feature_ls]
       
    if df_lable is not None:
        sample_name = np.array(genes_samples_mat.index)
        sample_lable = df_lable.loc[sample_name]
        sample_name = np.array([ sample_name[i]+' (%s)'%round(sample_lable[i],2) for i in range(len(sample_name)) ])
        genes_samples_mat.index = sample_name

    %R library(pheatmap)
    %R -i genes_samples_mat
    %R p_res = pheatmap(genes_samples_mat)


###################
## wrap function ##
###################
def data_filter(which, excluded_ls=None, included_ls=None):
    '''
    Input
        which: filter object is gene or sample
        excluded_ls: list or np.array 正则表达式??
        included_ls: list of np.array

    '''

def import_data(genomic_alter_path, expr_mat_path, drug_response_path):
    drug_response_df = pd.read_csv(drug_response_path,index_col=0)
    expr_mat = pd.read_csv(expr_mat_path,index_col=0)

    # remove deficiency data
    sample_list = list(set(drug_response_df.index) & set(expr_mat.columns))
    gene_list = list(expr_mat.columns)

    expr_mat = expr_mat.T.loc[sample_list]
    drug_response_df = drug_response_df.loc[sample_list]

    return {'expr_mat':expr_mat,
            'drug_response_df':drug_response_df,
            'sample_list':sample_list,
            'gene_list':gene_list
        }


def do_bootstrap_en(expr_mat, drug_response, resample_size, 
                    n_resample, replace_resample_flag, log_flag, 
                    model_max_iter, en_score_threshold, verbose=1):
    '''
    # input:
        resample_size: resample的结果的样本数
        n_resample: resample的次数,也就是获得多少组resample样本
        replace_resample_flag: 1/0 放回还是不放回抽样

        feature_number_threshold : 规定有多少feature,这次计算才算有效 [这里不用]
        iter_max : [这里不用]

        model_max_iter = 5000 : en模型最大的迭代数

        expr_mat: DataFrame, columns:gene, index:cell line
        drug_response: Series, index:cell line  不带药名 !!! [待解决]
        log_flag: 1/0 do np.log(x+1)

        en_score_threshold: paper said it shuold be 0.6

        verbose: 1/0 print verbose information

    # output:
        return {'model_ls':model_ls,
            'feature_dict':feature_dict,
            'feature_df':feature_df,
            'selected_feature_df':selected_feature_df,
            'resample_list':resample_list }   

    '''
    #############################
    ## do bootstrap resampling ##
    #############################
    #resample_size = 38
    #n_sample = 38
    #n_resample = 10
    n_sample = len(expr_mat)

    resample_list = []

    %R -i n_sample
    %R -i resample_size
    %R -i replace_resample_flag
    for i in range(n_resample):
        
        rs_result = %R sample(n_sample,size=resample_size,replace=replace_resample_flag)
        resample_list.append(rs_result)
        
    resample_list = np.array(resample_list)-1
    #############
    ## do iter ##
    #############
    if verbose:
        print('#Start# fit en model with use of %s groups of bootstrap samples...'%n_resample)

    # model_ls and feature_matrix
    model_ls = []
    feature_dict = {}

    for i_resample in range(len(resample_list)):
        if verbose:
            print('doing #%s iteration...'%(i_resample+1))
        ######################
        # get and split data #
        ######################
        gene_list = np.array(expr_mat.columns)

        if log_flag == 1:
            X_train = np.log(np.array( expr_mat.iloc[resample_list[i_resample]] )+1)
        else:
            X_train = np.array( expr_mat.iloc[resample_list[i_resample]] )
        
        y_train = np.array( drug_response.iloc[resample_list[i_resample]] )

        
        #########
        # model #
        #########
        en_res = do_elastic_net(X=X_train, y=y_train, 
                                all_feature_ls=gene_list, alphas=None, cv=10, 
                                random_state=random.randint(1,100), selection='random', 
                                max_iter=model_max_iter, verbose=verbose)
        model_ls.append(en_res)
        
        ################
        # feature dict #
        ################
        for f in en_res['feature_ls']:
            if f in feature_dict:
                feature_dict[f].append(en_res['feature_weight_df'].loc[f,'Weight'])
            else:
                feature_dict[f] = [en_res['feature_weight_df'].loc[f,'Weight']]
    if verbose:
        print('#End# totally find %s potential features'%len(feature_dict)) 

    #################################################################
    ## select feature with EN score larger than en_score_threshold ##
    #################################################################
    ## arrange feature from bootstrap_en ##
    for key in feature_dict.keys():
        # calculate
        n_posi = sum(np.array(feature_dict[key]) > 0)
        n_nega = sum(np.array(feature_dict[key]) < 0)        
        feature_dict[key] = ([np.array(feature_dict[key]), n_posi, n_nega])    
    
    feature_df = pd.DataFrame(feature_dict).T
    feature_df.columns = ['coef list','positive','negative']

    ## select feature according to en score ##
    selected_feature_ls = []
    for f in feature_df.index:
        en_score = max(feature_df.loc[f,'positive'], feature_df.loc[f,'negative'])/n_resample
        # print('feature is %s , the EN_score is %s'%(f, en_score))
        if en_score >= en_score_threshold:
            selected_feature_ls.append((f,en_score))
    selected_feature_df = pd.DataFrame(selected_feature_ls, columns=['feature','en score']).set_index(['feature'])

    print('Totally find %s feature with use of bootstrap+elastic net method !!\n===\nEND\n===')

    return {'model_ls':model_ls,
            'feature_dict':feature_dict,
            'feature_df':feature_df,
            'selected_feature_df':selected_feature_df,
            'resample_list':resample_list }


#################
# maybe useless #
def do_resample_en(expr_mat, drug_response, test_size, 
                    n_resample, replace_resample_flag, log_flag, 
                    model_max_iter, en_score_threshold, verbose=1):
    '''
    # input:
        expr_mat: DataFrame , columns is gene , index is sample
        drug_response: Series
        log_flag: 1/0 do np.log(x+1) default:1

        test_size: the percentage of test set
        iter_time: 模型迭代多少遍, 类似于do_bootstrap_en()中的 n_resample
        : 

    '''
    feature_number_threshold = 20
    iter_time = 100


    # model_ls and feature_matrix
    # resample_ls = []
    model_result_ls = []
    feature_dict = {}

    for i in range(iter_max):
        print('#%s iter time...'%(i+1))
        ######################
        # get and split data #
        ######################
        gene_list = np.array(expr_mat.columns)
        divide_res = train_test_set(df_samlple=expr_mat, df_lable=drug_response, test_size=0.2, return_type='data')

        if log_flag is 1:
            X_train = np.log( np.array(divide_res['df_samlple_train'])+1 )
            X_test = np.log( np.array(divide_res['df_samlple_test'])+1 )
        else:
            X_train = np.array(divide_res['df_samlple_train'])
            X_test = np.array(divide_res['df_samlple_test'])

        y_train = np.array(divide_res['df_label_train'])
        y_test = np.array(divide_res['df_label_test'])

        
        #########
        # model #
        #########
        en_res = do_elastic_net(X=X_train, y=y_train, all_feature_ls=gene_list, alphas=None, cv=10, 
                                random_state=random.randint(1,100), selection='random', max_iter=5000 )
        model_result_ls.append(en_res)

        ################
        # feature dict #
        ################
        for f in en_res_2['feature_ls']:
            if f in feature_dict:
                feature_dict[f].append(en_res['feature_weight_df'].loc[f,'Weight'])
            else:
                feature_dict[f] = [en_res['feature_weight_df'].loc[f,'Weight']]

    return {'feature_dict':feature_dict, 
            'model_result_ls':model_ls,
            'resample_list':resample_list}    
    
#################

def do_simple_en(expr_mat, drug_response, log_flag, 
                en_max_iter, feature_number_threshold, verbose=1):
    '''
    expr_mat   row: gene  columns:sample
    drug_response  Series
    '''
    gene_list = np.array(expr_mat.index)


    if log_flag is 1:
        X_train = np.log( np.array(expr_mat.T)+1 )
    else:
        X_train = np.array(expr_mat.T)

    y_train = np.array(drug_response)

    
    #########
    # model #
    #########
    en_res = do_elastic_net(X=X_train, y=y_train, all_feature_ls=gene_list, alphas=None, cv=10, 
                            random_state=random.randint(1,100), selection='random', max_iter=en_max_iter )    

    ##########
    #
    ##########
    importance_feature = abs(en_res['feature_weight_df']).sort_values('Weight').index[-feature_number_threshold:]


    mat_for_plot = expr_mat.T.loc[:,importance_feature].iloc[140:190]

    %R library(pheatmap)
    %R -i mat_for_plot
    %R p_res = pheatmap(mat_for_plot)

    order = %R p_res$tree_row$order
    labels = %R p_res$tree_row$label
    plt.plot(all_drug_sen_df.loc[drug].set_index('lab_id').loc[drug_expr_sample_ls,'auc'].loc[labels[order-1]])



def do_simple_en_old(expr_mat, drug_response, log_flag, test_size, 
                en_max_iter, feature_number_threshold, verbose=1):
    '''
    Input:
        expr_mat: DataFrame, columns:gene, index:cell line
        drug_response: Series
        log_flag: 1/0

        test_size: percentage of test set
        iter_max
        3en_cv
        #en_random_state
        #en_selection
        en_max_iter
        feature_number_threshold

        verbose: 1/0 default is 1

    Output:

    '''
    
    # feature_number_threshold = 20
    # iter_max = 100


    # model_ls and feature_matrix
    model_ls = []
    feature_dict = {}

    #for i in range(iter_max):
    print('#%s iter time...'%(i+1))
    ######################
    # get and split data #
    ######################
    gene_list = np.array(expr_mat.columns)
    divide_res = train_test_set(df_samlple=expr_mat, df_lable=drug_response, 
                                test_size=test_size, return_type='data')

    if log_flag is 1:
        X_train = np.log( np.array(divide_res['df_samlple_train'])+1 )
        X_test = np.log( np.array(divide_res['df_samlple_test'])+1 )
    else:
        X_train = np.array(divide_res['df_samlple_train'])
        X_test = np.array(divide_res['df_samlple_test'])

    y_train = np.array(divide_res['df_label_train'])
    y_test = np.array(divide_res['df_label_test'])

    
    #########
    # model #
    #########
    en_res = do_elastic_net(X=X_train, y=y_train, all_feature_ls=gene_list, alphas=None, cv=10, 
                            random_state=random.randint(1,100), selection='random', max_iter=en_max_iter )
        #if len(en_res['feature_ls']) < feature_number_threshold:
        #    break

    '''
    ###########
    # heatmap #
    ###########
    if log_flag:
        tran = 'log'
    else:
        tran = None

    feature_pheatmap(expr_mat, en_res['feature_ls'], tran=tran, df_lable=drug_response)

    # only the test data heatmap
    feature_pheatmap(expr_mat.iloc[divide_res['test_index_ls']], en_res['feature_ls'], tran=tran, 
                        df_lable=drug_response)

    #########################
    # Do linear and compare #
    #########################
    predict_y_test = en_res['model'].predict(X_test)

    # get train/test mat with only lasso feature , lable use the same as above
    if log_flag is 1:
        X_train_for_linear = np.log(np.array(expr_mat.loc[divide_res['train_ls'],en_res['feature_ls']])+1)
        X_test_for_linear = np.log(np.array(expr_mat.loc[divide_res['test_ls'],en_res['feature_ls']])+1)
    else:
        X_train_for_linear = np.array(expr_mat.loc[divide_res['train_ls'],en_res['feature_ls']])
        X_test_for_linear = np.array(expr_mat.loc[divide_res['test_ls'],en_res['feature_ls']]) 

    reg_linear = LinearRegression().fit(X_train_for_linear, y_train)
    predict_y_test_linear = reg_linear.predict(X_test_for_linear)

    # plot compare graph
    argsort_res = np.argsort(y_test)

    plt.figure(figsize=(20,6))

    plt.plot([y_test[i] for i in argsort_res],'g-')
    plt.plot([predict_y_test[i] for i in argsort_res],'b^')
    plt.plot([predict_y_test_linear[i] for i in argsort_res],'ys')
    '''
    return {'model_result':en_res, 'divide_res':divide_res}




from sklearn.ensemble import RandomForestClassifier
def do_random_forest(X, y, all_feature_ls, max_depth, random_state, verbose=1):
    '''
    X
    y
    all_feature_ls

    max_depth
    random_state
    '''
    regr = RandomForestRegressor(max_depth=max_depth, random_state=random_state)
    regr.fit(X, y)

    predict_y = regr.predict(X)

    feature_ls = []
    feature_weight_df = []
    for i in range(len(regr.feature_importances_)):
        if regr.feature_importances_[i]!=0:
            feature_ls.append(all_feature_ls[i])
            feature_weight_df.append([all_feature_ls[i], regr.feature_importances_[i]])
    feature_weight_df = pd.DataFrame(feature_weight_df,columns=['Feature','Weight']).set_index(['Feature'])
    
    if verbose:  
        print('Random Forest keep %s features...'%len(feature_ls))
    
    return {'model':regr, 'feature_ls':feature_ls, 
            'feature_weight_df':feature_weight_df, 
            'X':X, 'y':y}

def do_simple_random_forest(expr_mat, drug_response, log_flag, test_size, iter_max, 
                max_depth, feature_number_threshold, verbose=1):
    '''
    feature_number_threshold： 根据这个threshold取importance拍前的feature来输出，而不是重新开始

    Input:
        expr_mat: DataFrame, columns:gene, index:cell line
        drug_response: Series
        log_flag: 1/0

        test_size: percentage of test set
        iter_max
        3en_cv
        #en_random_state
        #en_selection
        en_max_iter
        feature_number_threshold

        verbose: 1/0 default is 1

    Output:
    '''
    # feature_number_threshold = 20
    # iter_max = 100

    # model_ls and feature_matrix
    model_ls = []
    feature_dict = {}

    for i in range(iter_max):
        print('#%s iter time...'%(i+1))
        ######################
        # get and split data #
        ######################
        gene_list = np.array(expr_mat.columns)
        divide_res = train_test_set(df_samlple=expr_mat, df_lable=drug_response, 
                                    test_size=test_size, return_type='data')

        if log_flag is 1:
            X_train = np.log( np.array(divide_res['df_samlple_train'])+1 )
            X_test = np.log( np.array(divide_res['df_samlple_test'])+1 )
        else:
            X_train = np.array(divide_res['df_samlple_train'])
            X_test = np.array(divide_res['df_samlple_test'])

        y_train = np.array(divide_res['df_label_train'])
        y_test = np.array(divide_res['df_label_test'])

        #########
        # model #
        #########
        rf_res = do_random_forest(X=X_train, y=y_train, all_feature_ls=gene_list, 
                                max_depth=max_depth, random_state=random.randint(1,100))
        if len(rf_res['feature_ls']) < feature_number_threshold:
            break

    ###########
    # heatmap #
    ###########
    if log_flag:
        tran = 'log'
    else:
        tran = None

    feature_pheatmap(expr_mat, rf_res['feature_ls'], tran=tran, df_lable=drug_response)

    # only the test data heatmap
    feature_pheatmap(expr_mat.iloc[divide_res['test_index_ls']], rf_res['feature_ls'], tran=tran, 
                        df_lable=drug_response)

    #########################
    # Do linear and compare #
    #########################
    predict_y_test = rf_res['model'].predict(X_test)

    # get train/test mat with only lasso feature , lable use the same as above
    if log_flag is 1:
        X_train_for_linear = np.log(np.array(expr_mat.loc[divide_res['train_ls'],rf_res['feature_ls']])+1)
        X_test_for_linear = np.log(np.array(expr_mat.loc[divide_res['test_ls'],rf_res['feature_ls']])+1)
    else:
        X_train_for_linear = np.array(expr_mat.loc[divide_res['train_ls'],rf_res['feature_ls']])
        X_test_for_linear = np.array(expr_mat.loc[divide_res['test_ls'],rf_res['feature_ls']]) 

    reg_linear = LinearRegression().fit(X_train_for_linear, y_train)
    predict_y_test_linear = reg_linear.predict(X_test_for_linear)

    # plot compare graph
    argsort_res = np.argsort(y_test)

    plt.figure(figsize=(20,6))

    plt.plot([y_test[i] for i in argsort_res],'g-')
    plt.plot([predict_y_test[i] for i in argsort_res],'b^')
    plt.plot([predict_y_test_linear[i] for i in argsort_res],'ys')



    return {'model_result':rf_res}


#############
# new
#############

def train_model(expr_mat, drug_response, model, sample_ls=None):
    '''
    expr_mat: DataFrame, columns:sample, index:gene
    drug_response: Series
    log_flag: 1/0
    '''
    # 获取有用的样本
    usrful_sample_set = set(expr_mat.columns) & set(drug_response.index)
    if sample_ls != None:
        usrful_sample_set = usrful_sample_set & set(sample_ls)
    # 选择使用的模型，训练模型
    if model = 'en':
        do_res = do_simple_en(expr_mat, 
                    drug_response, 
                    log_flag, 
                    en_max_iter, 
                    feature_number_threshold, 
                    verbose=1)        




if __name__ == '__main__':
    genomic_alter_path = ''
    expr_mat_path = '/y/home/lgh/lgh/work/drug_sensitivity/data/all_breast_cell_com_geneexpression_norm.csv'
    drug_response_path = '/y/home/lgh/lgh/work/drug_sensitivity/data/ic50_formatted.csv'

    import_data_res = import_data(genomic_alter_path=None,
            expr_mat_path=expr_mat_path,drug_response_path=drug_response_path)


    expr_mat = import_data_res['expr_mat']
    drug_response_df = import_data_res['drug_response_df']
    

