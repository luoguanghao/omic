import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


tp53_func = '../data/tp53_iarc/geneVariationIARC TP53 Database, R20.txt'
tp53_func_df = pd.read_csv(tp53_func, sep='\t')

mut_df = pd.read_csv('../data/CCLE_mutations.csv')


mut_df_anno_fun = pd.merge(
    mut_df.loc[mut_df['Hugo_Symbol']=='TP53'],tp53_func_df, left_on='Protein_Change', right_on='ProtDescription')

mut_df_anno_fun.shape










