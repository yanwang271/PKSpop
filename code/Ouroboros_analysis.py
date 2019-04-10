#!/usr/bin/env python3

"""
Run Ouroboros analysis and select the result with best LLH
"""

import json
import os
import subprocess

def ouroboros_analysis(info_dict):
    '''
    '''
    query_c = info_dict['c_query_path']
    query_n = info_dict['n_query_path']
    if info_dict['Ouroboros_int_frac']:
        param_ls = info_dict['Ouroboros_int_frac']
    else:
        param_ls = [0.90, 0.80] #default 
    n_repeat = info_dict['n_repeat']
    ouro_path = info_dict['Ouroboros_path']
    
    ### Run Ouroboros repeatly with defined int_frac parameter
    llh_param = ouroboros_best_result(ouro_path, query_c, param_ls, n_repeat)
    
    ### Find the Ouroboros result with the best LLH
    param_idx, repeat_idx = best_llh(llh_param)
    param = param_ls[param_idx]
    
    info_dict.update({'best_param':param,
                      'best_repeat':repeat_idx})
                      
    return info_dict

def ouroboros_best_result(ouro_path, query_c, param_ls, n_repeat):
    '''
    Run Ouroboros repeatly with different int_frac parameter
    '''    
    llh_param = []
    for p in param_ls:
        llh_repeat = []
        for repeat in range(n_repeat):
            json_name = query_c.replace('.fasta',f'_{p}_soft_warm_{repeat}.json')\
                .replace('/c-dd','/Ouroboros')
            if not os.path.exists(json_name):
                create_json_file(json_name, query_c, p)
                ouro_res_dir = json_name.replace('.json','')   
            
            run_Ouroboros(json_name, ouro_path)
            llh_repeat.append(extract_llh(ouro_res_dir))
        llh_param.append(llh_repeat)
    return llh_param

def run_Ouroboros(json_inpt, ouro_path):
    '''
    Run Ouroboros analysis
    '''
    run_cmd = f'python {ouro_path}/code/src/run_analysis.py {json_inpt}'
    subprocess.check_call(run_cmd.split(' '))

def extract_llh(res_directory): 
    '''
    '''
    llh_fl = open(f'{res_directory}/all_total_llhs.csv')
    llh = llh_fl.readlines()[-1].strip()
    return float(llh)

def best_llh(llh_mtx):
    '''
    '''
    param_max = [max(ls) for ls in llh_mtx]
    param_idx = param_max.index(max(param_max))
    repeat_idx = llh_mtx[param_idx].index(max(param_max))
    return param_idx, repeat_idx

def create_json_file(json_name, query_fl, param):
    '''
    Create json file that containing Ouroboros input information
    '''
    json_fl = open(json_name,'w')
    inpt_info = {
                 "io": json_name.replace('.json',''),
                 "mode": "soft",
                 "init": "warm",
                 "msa1": query_fl.replace('/c-dd','/n-dd'),
                 "msa2": query_fl,
                 #"test": True,
                 #"int_limit": 222,
                 "int_frac": param,
                 "gap_threshold": 0.5,
                 "n_jobs": 20
                }
    json_fl.write(json.dumps(inpt_info, indent=2))
    json_fl.close()
