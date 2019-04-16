#!/usr/bin/env python3

"""
Using the interaction probabilities resulted from Ouroboros analysis to
predict the protein order
"""

from pandas import DataFrame as df
from pandas import Series
import numpy as np

def predict_order(info_dict):
    '''
    '''
    print(f'Predicting protein order...')
    output_dir = info_dict['output_path']
    old_groups = [[]]
    
    if info_dict['other_class_c'] and info_dict['other_class_n']:
        old_groups = add_notclass1_pairs(info_dict['other_class_c'],\
                            info_dict['other_class_n'])

    find_next_max = True
    res_fl = open(f'{output_dir}/prediction_result.txt','w')
    
    info_dict['prob_mtx'] = info_dict['prob_mtx'].astype('float64')
    prob_mtx = info_dict['prob_mtx'].replace(np.nan, -1)
    prob_mtx = remove_same_pair(prob_mtx)
    
    while find_next_max:
        max_val_pos = find_max_prob(prob_mtx)

        if len(max_val_pos) == 1:
            new_groups = add_pair_to_groups(old_groups, max_val_pos[0])
            prob_mtx = remove_val_mtx(prob_mtx, max_val_pos)
            old_groups = new_groups[:]
    
        else:
            con_pairs, uncon_pairs = grouping_pairs(max_val_pos)
            for p in uncon_pairs:
                new_groups = add_pair_to_groups(old_groups, p)
                old_groups = new_groups[:]
            prob_mtx = remove_val_mtx(prob_mtx, uncon_pairs)
            
            for p_p in con_pairs:
                new_groups = []
                for p in p_p:
                    new_subgroups = add_pair_to_groups(old_groups, p)
                    new_groups += new_subgroups
                old_groups = set_list_of_list(new_groups)
                prob_mtx = remove_val_mtx(prob_mtx, p_p)
        
        if prob_mtx.max().max() == -1:
            find_next_max = False
        
        res_info = df(columns=['All_proteins','Correct_start_protein','Correct_end_protein'])
        for group in old_groups:
            line = assemble_line(group)
            if not line:
                break
            
            ###overall_prob = cal_overall_prob(info_dict['prob_mtx'], line)
            order = ','.join(line)
            ###res_info.loc[(order,'Overall_prob')] = overall_prob
            if complete_line(line, info_dict['protein_id']):            
                
                res_info.append(Series(name=order))
                res_info.loc[(order,'All_proteins')] = True
                
                if info_dict['start_pro']: 
                    if correct_start(line, info_dict['start_pro']):
                        res_info.loc[(order, 'Correct_start_protein')] = True
                    else:
                        res_info.loc[(order, 'Correct_start_protein')] = False
                
                if info_dict['end_pro']:
                    if correct_end(line, info_dict['end_pro']):
                        res_info.loc[(order, 'Correct_end_protein')] = True
                    else:
                        res_info.loc[(order, 'Correct_end_protein')] = False
                
                if res_info.at[order,'Correct_end_protein'] != False and \
                   res_info.at[order,'Correct_start_protein'] != False:
                       find_next_max = False
                
            else:
                res_info.loc[(order,'All_proteins')] = False
                        
    res_info = res_info[res_info['Correct_start_protein'] != False]
    res_info = res_info[res_info['Correct_end_protein'] != False]
    
    complete_lines = res_info[res_info['All_proteins'] == True]
    rows = list(complete_lines.index)
    '''
    if len(rows) >= 1:
        max_overall_prob = complete_lines['Overall_prob'].max()
        for r in rows:
            if complete_lines.at[r,'Overall_prob'] == max_overall_prob:
                pred_order = r
    else:
        max_overall_prob = res_info['Overall_prob'].max()
        for r in list(res_info.index):
            if res_info.at[r,'Overall_prob'] == max_overall_prob:
                pred_order = r
    '''
    if len(rows) >= 1:
        res_fl.write('The predicted order is:\n')
        for pred_order in rows:
            print(f'The predicted order is {pred_order}')
            res_fl.write(pred_order+'\n')
    else:
        res_fl.write('No prediction')
        print('No prediction')
    res_fl.close()


def add_notclass1_pairs(not_c, not_n):
    not_pairs = []
    for c in not_c:
        for n in not_n:
            not_pairs.append((c,n))
    
    groups = []
    for i in range(len(not_c)):
        p1 = not_pairs[i]
        group = [p1,]
        for j in range(i+1,len(not_pairs)):
            p2 = not_pairs[j]
            if not wether_conflict([p1],p2):
                group.append(p2)
        groups.append(group)
    return groups
                
def remove_same_pair(mtx):
    '''
    Remove the interaction probabilities of pairs from the same protein
    '''
    columns = list(mtx)
    rows = list(mtx.index)
    for c in columns:
        for r in rows:
            if c == r:
                mtx = remove_val_mtx(mtx, [(c,c)])
    return mtx

def find_max_prob(prob_mtx):
    '''
    Find the highest interaction probability in the matrix
    '''
    max_val = prob_mtx.max().max()
    pos_val = idx_val(prob_mtx, max_val)
    return pos_val

def idx_val(mtx, val):
    '''
    Find the positions of a value in the matrix
    '''
    pos_val = []
    columns = list(mtx)
    rows = list(mtx.index)
    for c in columns:
        for r in rows:
            if c != r and mtx.at[r, c] == val:
                pos_val.append((r,c))
    return pos_val

def add_pair_to_groups(groups, pair):
    '''
    Add interacting protein pair to existing pair groups that are not 
    conflict with the pair
    
    Arguements
    groups:     list, consists of group lists 
                (group: list, consist of interacting protein pair tuples)
    pair:       tuple, interacting protein pair (c_docking_domain, n_docking_domain)
    '''
    new_groups = []
    for group in groups:
        if len(group) < 1:
            new_group = group + [pair]
        else:
            conflict = wether_conflict(group, pair)
            if not conflict:
                new_group = group + [pair]
            else:
                new_group = group
        new_groups.append(new_group)
    return new_groups
    
def grouping_pairs(pair_ls):

    '''
    For protein pairs with the same interacting prob, find the devide 
    the conflict pairs and unconflict pairs
    '''
    conflict_pair = []
    unconflict_pair = []
    for p1 in pair_ls:
        conflict = False
        for p2 in pair_ls:
            if p1 != p2:
                if wether_conflict([p2],p1):
                    conflict = True
                    if (p2,p1) not in conflict_pair:
                        conflict_pair.append((p1,p2))
        if not conflict:
            unconflict_pair.append(p1)
    return conflict_pair, unconflict_pair
            
def remove_val_mtx(mtx, position):
    '''
    '''
    for p in position:
        r, c = p
        mtx.at[r, c] = -1
    return mtx
                
def set_list_of_list(ls):
    '''
    '''
    res_ls = []
    for item in ls:
        if item not in res_ls:
            res_ls.append(item)
    return res_ls
    
def wether_conflict(pair_ls, pair):
    '''
    Judge wether the pair is conflict with the pairs in pair_ls
    '''
    c_ls = []
    n_ls = []
    c, n = pair
    for p in pair_ls:
        c_ls.append(p[0])
        n_ls.append(p[-1])
        if p == (n,c):
            return True
    
    if c in c_ls:
        return True
    elif n in n_ls:
        return True
    else:
        pair_ls_ = pair_ls[:]
        loop = wether_loop(pair_ls_, pair)
        return loop

def wether_loop(pair_ls, pair):
    c, n = pair
    c_ls = [i[0] for i in pair_ls]
    n_ls = [i[1] for i in pair_ls]
    length = len(c_ls)
    
    if c in c_ls and n in n_ls:
        for i in range(length):
            if c in c_ls[i]:
                new_c = n_ls[i]
                break
        for j in range(length):
            if n in n_ls[j]:
                new_n = c_ls[j]
                break
        return next_check(pair_ls, new_c, new_n, i, j)        
    
    if c in n_ls and n in c_ls:
        for i in range(length):
            if c in n_ls[i]:
                new_n = c_ls[i]
                break
        for j in range(length):
            if n in c_ls[j]:
                new_c = n_ls[j]
                break
        return next_check(pair_ls, new_c, new_n, i, j)        
    else:
        return False
        
def next_check(pair_ls, new_c, new_n, i, j):
        new_pair = (new_c, new_n)
        if new_pair in pair_ls:
            return True
        elif new_c == new_n:
            return True
        else:
            rm1 = pair_ls[i]
            rm2 = pair_ls[j]
            pair_ls.remove(rm1)
            pair_ls.remove(rm2)
            loop = wether_loop(pair_ls, new_pair)
            return loop        

def assemble_line(pair_ls):
    '''
    Assemble the line from the predicted interacting protein pairs
    '''
    line = []
    c_ls = []
    n_ls = []
    for p in pair_ls:
        c_ls.append(p[0])
        n_ls.append(p[-1])
    
    loop = True          
    for c in c_ls:
        if c not in n_ls:
            line.append(c)
            loop = False
            break
    if loop: 
        return False        
    
    while c:
        try:
            c_idx = c_ls.index(c)
            n = n_ls[c_idx]
            line.append(n)
            c = n
        except ValueError:
            c = False
    return tuple(line)
        
def complete_line(line, protein_ls):
    '''
    Judge wether the assembled line is complete
    '''
    if len(line) == len(protein_ls):
        return True
    else:
        return False
    
def correct_start(line, start_pro):
    '''
    Judge wether the assembled line start/end of the correct protein 
    '''
    if line[0] == start_pro:
        return True
    else:
        return False
    
def correct_end(line, end_pro):
    '''
    '''
    if line[-1] == end_pro:
        return True
    else:
        return False

def cal_overall_prob(prob_mtx, line):
    '''
    Calculate the the oveall probability by adding the the probability
    of all the pairs
    '''
    prob_mtx = prob_mtx.replace(np.nan, 0)
    overall_prob = 0
    for i in range(len(line)-1):
        c = line[i]
        n = line[i+1]
        prob = prob_mtx.at[c,n]
        overall_prob = overall_prob + prob
    return overall_prob
    
