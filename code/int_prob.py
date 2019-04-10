#!/usr/bin/env python3

"""
Generate pairwise interaction probability matrix
"""

from pandas import DataFrame as df
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import subprocess
import os
import json
import csv

def prob_mtx(info_dict):
    '''
    '''    
    output_dir = info_dict['output_path']
    pro_ls = info_dict['protein_id']
    query_c = info_dict['c_query_path']
    query_n = info_dict['n_query_path']
    param = info_dict['best_param']
    repeat = info_dict['best_repeat']
    end_pro = info_dict['end_pro']
    start_pro = info_dict['start_pro']
    
    probs = extract_prob(query_c, param, repeat)
    c_pro = pro_inpt_order(query_c)
    n_pro = pro_inpt_order(query_n)
    
    prob_mtx = create_prob_mtx(probs, c_pro, n_pro, pro_ls, output_dir)
    prob_mtx_heatmap(prob_mtx, output_dir)
    
    other_class_c = find_other_class(pro_ls, c_pro, end_pro)
    other_class_n = find_other_class(pro_ls, n_pro, start_pro)
    print(other_class_c, other_class_n)
    
    
    info_dict.update({'other_class_c':other_class_c,
                      'other_class_n':other_class_n,
                      'prob_mtx':prob_mtx})

    return info_dict

def extract_prob(query_c, param, repeat):
    '''
    '''
    ouro_res_dir = query_c.replace('.fasta',f'_{param}_soft_warm_{repeat}')\
                .replace('/c-dd','/Ouroboros')
    all_labels = list(csv.reader(open(f'{ouro_res_dir}/labels_per_iter.csv',\
                 newline = ''), delimiter = ' '))
    labels = all_labels[523:]
    probs = [format(float(i[-1])) for i in labels]
    return probs
    
def pro_inpt_order(query):
    '''
    '''
    record = list(SeqIO.parse(query,'fasta'))[523:]
    pro_order = []
    for r in record:
        pro_order.append(r.name)
    return pro_order
    
def create_prob_mtx(probs, c_pro, n_pro, pro_ls, output_dir):
    '''
    Organize the interaction probability from Ouroboros analysis to a
    dataframe
    '''
    ### The matrix contain all the proteins in the assembly line
    print(f'Creating interaction probability matrix...')
    prob_mtx = df(index = pro_ls, columns = pro_ls)
    
    for i in range(len(probs)):
        c_term = c_pro[i]
        n_term = n_pro[i]
        prob_mtx.loc[c_term][n_term] = probs[i]
    prob_mtx.to_csv(path_or_buf=f'{output_dir}/int_prob_mtx.csv')
    return prob_mtx

def prob_mtx_heatmap(prob_mtx, output_dir):
    '''
    Generate interaction probability heatmap
    '''
    prob_mtx = prob_mtx.astype(np.float64).apply(np.log10)
    prob_mtx = prob_mtx.applymap(minus)
    sns.set()
    
    fig, ax = plt.subplots(figsize=(5,5))
    ax = sns.heatmap(prob_mtx, cmap="YlGnBu", annot=True,\
                     annot_kws={'size':10}, ax=ax) #
    fig.savefig(f'{output_dir}/int_prob_mtx.png')
    
def minus(val):
    '''
    '''
    return 0-val
    
def find_other_class(all_pro, class1_pro, start_end):
    '''
    Find the proteins that not belong to class 1
    '''
    other_class = []
    for p in all_pro:
        if p not in class1_pro and p not in start_end:
            other_class.append(p)
    return other_class
            
