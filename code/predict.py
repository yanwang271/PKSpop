#!/usr/bin/env python3

"""
Using the interaction probabilities resulted from Ouroboros analysis to
predict the protein order
"""

import pandas as pd
import numpy as np
from scipy import optimize

def predict_order(info_dict):
    '''
    Predict protein order from the interaction probability matrix and
    information of compatibility class 
    '''
    print(f'Predicting protein order...')
    output_dir = info_dict['output_path']
    pro_ls = info_dict['protein_id']
    
    res_order = open(f'{output_dir}/protein_order_prediction.txt','w')
    res_ppi = open(f'{output_dir}/protein_interaction_prediction.txt','w')

    info_dict['prob_mtx'] = info_dict['prob_mtx'].astype('float64')
    prob_mtx = info_dict['prob_mtx'].values
        
    prob_mtx[np.isnan(prob_mtx)] = 0
    prob_mtx = remove_same_pair(prob_mtx) 
    
    if info_dict['other_class_c'] and info_dict['other_class_n']:
       
        other_class_c = [pro_ls.index(p) for p in info_dict['other_class_c']]
        other_class_n = [pro_ls.index(p) for p in info_dict['other_class_n']]
        notclass1_pairs = add_notclass1_pairs(other_class_c, other_class_n)
        prob_mtx = remove_val_mtx(prob_mtx, notclass1_pairs, 1)
    
    prob_mtx = -prob_mtx

    row_ind, col_ind = optimize.linear_sum_assignment(prob_mtx)
    end = info_dict['end_pro']
    start = info_dict['start_pro']
    
    c_dd = [pro_ls[i] for i in list(row_ind)]
    n_dd = [pro_ls[i] for i in list(col_ind)]
    
    groups = list(zip(c_dd, n_dd)) 
    
    remove = []
    for pair in groups:
        if pair[1] == start:
            remove.append(pair)
        
        if pair[0] == end:
            remove.append(pair)
                
        if pair[0] == pair[1]:
            remove.append(pair)
    
    groups = list(set(groups) - set(remove))
    pairs = ['-'.join(p) for p in groups]
    res_ppi.write(f'The predicted interacting protein pairs (Cdd-Ndd) are:\
    \n{", ".join(pairs)}\n')
    
    line = assemble_line(groups)
    order = ','.join(line)
    print(f'The predicted order is:\n{order}')
    res_order.write(f'The predicted order is:\n{order}')
 
    res_order.close()
    res_ppi.close()


def add_notclass1_pairs(not_c, not_n):
    '''
    Paring DDs that do not belong to class I
    '''
    not_pairs = []
    for c in not_c:
        for n in not_n:
            not_pairs.append((c,n))
    return not_pairs
                
def remove_same_pair(mtx):
    '''
    Assign interaction probability = 0 to pairs from the same protein
    '''
    columns = mtx.shape[1]
    for c in range(columns):
        mtx = remove_val_mtx(mtx, [(c,c)], 0)
    return mtx
            
def remove_val_mtx(mtx, position, to_value):
    '''
    Assign a value to the given position in an matrix
    '''
    for p in position:
        r, c = p
        mtx[r, c] = to_value
    return mtx
                
def assemble_line(pair_ls):
    '''
    Assemble the line from the predicted interacting protein pairs
    '''
    line = []
    c_ls = [p[0] for p in pair_ls]
    n_ls = [p[-1] for p in pair_ls]
    '''
    c_ls = []
    n_ls = []
    
    for p in pair_ls:
        c_ls.append(p[0])
        n_ls.append(p[-1])
    '''
    
    loop = True          
    for c in c_ls:
        if c not in n_ls:
            line.append(c)
            loop = False
            break
    if loop: 
        return 'loop'        
    
    while c:
        try:
            c_idx = c_ls.index(c)
            n = n_ls[c_idx]
            line.append(n)
            c = n
        except ValueError:
            c = False
    return tuple(line)
