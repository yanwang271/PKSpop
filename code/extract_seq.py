#!/usr/bin/env python3

"""
Extract docking domain sequences from gbk file
Filter sequences and find the start or end domain
"""

from Bio import SeqIO
import warnings

def extract_seq(info_dict):
    '''
    Extract dd sequences from gbk file
    
    Arguements
    gbk_path:   str, the file path of gbk file of query gene cluster
    pro_ls:     list, the list of query protein identifiers
    id_category:   str, the category that the protein identifiers belong
    res_dir:    str, the path of directory that store data and output
    
    Return
    info_dict:  dict, dictionary that contain path to sequence fasta file,
                      start/end protein, query protein list
    '''
    gbk_path = info_dict['gbk_path']
    pro_ls = info_dict['protein_id']
    id_category = info_dict['id_category']
    res_dir = info_dict['result_path']

    #info_dict.update({'result_path':res_dir})
    c_seq_path = f'{res_dir}/c-dd_raw.fasta'
    n_seq_path = f'{res_dir}/n-dd_raw.fasta'
    c_dd = open(c_seq_path,'w')
    n_dd = open(n_seq_path,'w')
    info_dict.update({'c_seq_path':c_seq_path,
                      'n_seq_path':n_seq_path,
                      'start_pro':[],
                      'end_pro':[]})
    
    gbk_fl = open(gbk_path)
    record = list(SeqIO.parse(gbk_fl, 'genbank'))[0]
    
    tool = 'antismash'
    for f in record.features:
        if f.type == 'subregion':
            if 'mibig' in f.qualifiers['aStool']:
                tool = 'mibig'
                break
    
    for pro in pro_ls:
        print(f'Extracting sequences from {pro}...')
        
        if tool == 'antismash':
            c_seq, n_seq = find_dd(pro, record, id_category)
        elif tool == 'mibig':
            c_seq, n_seq = find_dd_mibig(pro, record, id_category)

        if c_seq == 'end' and n_seq == 'start':
            warnings.warn(f'{pro} do not have c and n term docking domain, \
            cannot place it into the assembly line.')
            info_dict['protein_id'].remove(pro)
            continue
            
        if c_seq == 'end':
            info_dict['end_pro'].append(pro)
        else:
            c_dd.write(f'>{pro}\n{c_seq}\n')
        
        if n_seq == 'start':
            info_dict['start_pro'].append(pro)
        else:
            n_dd.write(f'>{pro}\n{n_seq}\n')
    
    gbk_fl.close()
    c_dd.close()
    n_dd.close()
    
        
def find_dd(pro, record, id_category):
    '''
    Find n and c-term docking domain sequences on a protein
    '''
    n_seq = 'start'
    c_seq = 'end'
    for f in record.features:
        if f.type == 'CDS':
            try:
                f.qualifiers[id_category]
            except KeyError:
                continue
            if pro in f.qualifiers[id_category]:
                sec_met = f.qualifiers['sec_met']
                sec_met = [i for i in sec_met if 'ASF-prediction' not in i]
                for sec in sec_met:
                    if sec.startswith('NRPS/PKS Domain:'):
                        sec_info = sec.split(' ')
                        
                        ### Find n_term docking domain
                        if 'PKS_KS' in sec_info[2]:
                            last_sec = sec_met[sec_met.index(sec)-1]
                            if 'NRPS/PKS subtype' in last_sec or \
                               'Docking' in last_sec:
                                ks_start, ks_end = parse_loc(sec_info[3])
                                #if ks_start >= 26:
                                if ks_start > 0:
                                    n_seq = f.qualifiers['translation']\
                                            [0][:ks_start]
                        ### Find c_term docking domain
                        if 'ACP' in sec_info[2]:
                            try:
                                next_sec = sec_met[sec_met.index(sec)+1]
                                if 'Docking' in next_sec:
                                    acp_start, acp_end = parse_loc(sec_info[3])
                                    c_seq = f.qualifiers['translation']\
                                            [0][acp_end:]
                                elif 'Thioesterase' in next_sec:
                                    c_seq = 'end'
                                else:
                                    continue
                            except IndexError:
                                acp_start, acp_end = parse_loc(sec_info[3])
                                c_seq = f.qualifiers['translation']\
                                        [0][acp_end:]
        
                return c_seq, n_seq

def find_dd_mibig(pro, record, id_category):
    '''
    Find n and c-term docking domain sequences on a protein in mibig gbk file
    '''
    n_seq = 'start'
    c_seq = 'end'
    for f in record.features:
        if f.type == 'CDS':
            try:
                f.qualifiers[id_category]
            except KeyError:
                continue
            if pro in f.qualifiers[id_category]:
                sec_met = f.qualifiers['NRPS_PKS']
                sec_met = [i for i in sec_met if 'Domain' in i]
                for sec in sec_met:
                    sec_info = sec.split(' ')
                    
                    ### Find n_term docking domain
                    if 'PKS_KS' in sec_info[1]:
                        if sec_met.index(sec) == 0:
                            ks_start, ks_end = parse_loc(sec_info[2])
                            #if ks_start >= 26:
                            if ks_start > 0:
                                n_seq = f.qualifiers['translation']\
                                        [0][:ks_start]
                        else:
                            last_sec = sec_met[sec_met.index(sec)-1]
                            if 'Docking' in last_sec:
                                ks_start, ks_end = parse_loc(sec_info[2])
                                #if ks_start >= 26:
                                if ks_start > 0:
                                    n_seq = f.qualifiers['translation']\
                                            [0][:ks_start]
                            else:
                                continue                            
                               
                    ### Find c_term docking domain
                    if 'ACP' in sec_info[1]:
                        try:
                            next_sec = sec_met[sec_met.index(sec)+1]
                            if 'Docking' in next_sec:
                                acp_start, acp_end = parse_loc(sec_info[2])
                                c_seq = f.qualifiers['translation']\
                                        [0][acp_end:]
                            elif 'Thioesterase' in next_sec:
                                c_seq = 'end'
                            else:
                                continue
                        except IndexError:
                            acp_start, acp_end = parse_loc(sec_info[2])
                            c_seq = f.qualifiers['translation']\
                                    [0][acp_end:]
        
                return c_seq, n_seq
    
def parse_loc(loc_info):
    '''
    '''
    start,end = list(map(int,\
                loc_info.strip('.').strip('()').split('-')))
    start = start
    end = end
    return start,end
    
    
def remove_invalid_pro(info_dict):
    '''
    Remove the protein if it do not have c and n termini
    '''
    for pro in info_dict['start_pro']:
        if pro in info_dict['end_pro']:
            warnings.warn(f'{pro} do not have c and n term docking \
            domain, cannot place it into the assembly line.')
            info_dict['protein_id'].remove(pro)
            info_dict['end_pro'].remove(pro)
            info_dict['start_pro'].remove(pro)

