#!/usr/bin/env python3

"""
Cluster and align the docking domain sequences
"""

import os
import subprocess 
from Bio import SeqIO
from extract_seq import remove_invalid_pro

def clustering(info_dict):
    '''
    Cluster the query sequences to 3 classes
    '''
    c_seq = info_dict['c_seq_path']
    n_seq = info_dict['n_seq_path']
    pro_ls = info_dict['protein_id']
        
    print(f'Clustering sequence...')
    c_class_seq = cluster_seq(c_seq,'PKSpop/data/hmm_profile/c_groups_hmmpf')
    n_class_seq = cluster_seq(n_seq,'PKSpop/data/hmm_profile/n_groups_hmmpf')
    
    for pro in pro_ls:
        if pro not in c_class_seq:
            if not info_dict['end_pro']:
                info_dict["end_pro"].append(pro)
                
        if pro not in n_class_seq:
            if pro not in info_dict['start_pro']:
                info_dict["start_pro"].append(pro)
    
    ### Check if there are proteins that do not have c and n termini
    remove_invalid_pro(info_dict)
    ### Check if there are more than one start/end protein
    wether_predict(info_dict)

    c_class1_fl = c_seq.replace('raw.fasta',\
                                 f'class_1.fasta')
    n_class1_fl = n_seq.replace('raw.fasta',\
                                 f'class_1.fasta')
    
    if not os.path.exists(c_class1_fl):
        raise Exception('There is no c-term class 1 docking domain in\
        this assembly line, cannot perform prediction')
    if not os.path.exists(n_class1_fl):
        raise Exception('There is no n-term class 1 docking domain in\
        this assembly line, cannot perform prediction')

    
    with open(c_class1_fl) as c_class1_seq:
        num_seq_c = len(list(SeqIO.parse(c_class1_seq,'fasta')))
    with open(n_class1_fl) as n_class1_seq:
        num_seq_n = len(list(SeqIO.parse(n_class1_seq,'fasta')))
    if num_seq_c != num_seq_c:
        raise Exception('Unequal number of class 1 c and n docking domain,\
                        cannot perform the protein order prediction')
                        
    return info_dict
    
    
def cluster_seq(query_fl, hmmpf_db):
    '''
    Cluster the query sequences to 3 classes
    Generate sequence fasta files of of 3 classes 
    '''
    hmmscan_out_path = query_fl.replace('raw.fasta',\
                                        'hmmscan_oupt.txt')
    ### Run hmmscan to asign each seq to its best match class
    run_hmmscan(hmmscan_out_path, hmmpf_db, query_fl)
    group, class_seq = parse_hmmscan_out(hmmscan_out_path)
    group_seq_fasta(group, query_fl)
    
    return class_seq
    

def msa(info_dict):
    '''
    Aligning sequences by HMMER
    '''
    c_seq = info_dict['c_seq_path']
    n_seq = info_dict['n_seq_path']

    ### Use class1 sequences to perform MSA
    c_inpt = c_seq.replace('raw.fasta',\
                              f'class_1.fasta')
    n_inpt = n_seq.replace('raw.fasta',\
                              f'class_1.fasta')
    ### Add a already cut sequence to help locate the cutting position                         
    add_position_helper(c_inpt, 'c')
    add_position_helper(n_inpt, 'n')
    ### Align sequeces against hmm profile                            
    c_align = c_inpt.replace('.fasta','_aln.afa')
    n_align = n_inpt.replace('.fasta','_aln.afa')
    print(f'Aligning sequences...')
    run_hmmalign(c_align, 'PKSpop/data/hmm_profile/c_group_1.hmm', c_inpt)
    run_hmmalign(n_align, 'PKSpop/data/hmm_profile/n_group_1.hmm', n_inpt)
    ### Identify the cutting position and cut the sequence
    c_start, c_end = find_cut_position(c_align)
    c_cut_fl = cut_seq(c_align, c_start, c_end)
    
    n_start, n_end = find_cut_position(n_align)
    n_cut_fl = cut_seq(n_align, n_start, n_end)
    
    info_dict.update({'c_align_path':c_cut_fl, 'n_align_path':n_cut_fl})
    return info_dict

def cut_seq(aln_fl, start_pos, end_pos):
    '''
    Cut the aligned sequences to obtain the part that used to perform 
    coevolutionary analysis
    '''
    cut_fl = aln_fl.replace('_aln.afa','.afa')
    cut_seq = open(cut_fl, 'w')
    aln_seq = open(aln_fl)
    record = list(SeqIO.parse(aln_seq, 'fasta'))
    for r in record:
        if 'positioning' not in r.name:
            seq = str(r.seq)[start_pos:end_pos]
            cut_seq.write(f'>{r.name}\n{seq}\n')
    aln_seq.close()
    cut_seq.close()
    return cut_fl
    
def add_position_helper(seq_fl, c_n):
    '''
    Add a already cut sequence to help locate the cutting position
    '''
    pos_seq = open('PKSpop/data/hmm_profile/positioning_helper.fasta')
    record = list(SeqIO.parse(pos_seq, 'fasta'))
    with open(seq_fl,'a') as seq:
        if c_n == 'c':
            r = record[0]
        else:
            r = record[1]
        seq.write(f'>{r.name}\n{str(r.seq)}\n')
    pos_seq.close()
    
def find_cut_position(aln_fl):
    '''
    Find the cutting position on the aligned sequences
    '''
    with open(aln_fl) as aln_seq:
        record = list(SeqIO.parse(aln_seq, 'fasta'))
        pos_seq = str(record[-1].seq)
        for i in range(len(pos_seq)):
            if pos_seq[i].isalpha() and not_alpha_s(pos_seq[:i]):
                start_pos = i
            elif pos_seq[i].isalpha() and not_alpha_s(pos_seq[i+1:]):
                end_pos = i+1
    return start_pos, end_pos
            
def not_alpha_s(string):
    '''
    Find wether there is a alpha in a string
    '''
    res = True
    for s in string:
        if s.isalpha():
            res = False
    return res
                
def run_hmmalign(aligned_fl, hmmfile, input_fl):
    '''
    Call HMMER
    multiple sequence alignment against the profile
    '''
    cmd = f'hmmalign -o {aligned_fl} --outformat afa {hmmfile} {input_fl}'
    subprocess.check_call(cmd.split(' '))

def run_hmmscan(hmmscan_out,hmmfile_db,query_fl):
    '''
    Search query sequences against hmmfile database, 
    to find out which family they belong
    '''
    cmd = f'hmmscan --tblout {hmmscan_out} {hmmfile_db} {query_fl}'
    subprocess.check_call(cmd.split(' '))
    
def parse_hmmscan_out(hmmscan_out):
    '''
    Find which group each sequence belong to and
    Store the group information into a dictionary: {group_name:[sequence_id],}
    '''
    scan = open(hmmscan_out)
    group = {}
    class_seq = []
    last_seq = ''
    for line in scan:
        if not line.startswith('#'):
            info = line.strip().split(' - ')
            info = [info[0].strip(),info[1].strip()]
            class_seq.append(info[1])
            if not info[0] in group:
                group.update({info[0]:[]})
            if info[1] != last_seq:
                group[info[0]].append(info[1])
                last_seq = info[1]
            else:
                last_seq = info[1]
    scan.close()
    return group, class_seq
    
def group_seq_fasta(group, query_fl):
    '''
    Write sequences from a group, which is detected by hmmscan, 
    to a new fasta file
    '''
    query_seq = open(query_fl)
    record = list(SeqIO.parse(query_seq,'fasta'))
    for name in group.keys():
        out_fl = query_fl.replace('raw.fasta',\
                                  f'class_{name[-1]}.fasta')
        fasta = open(out_fl,'w')
        for item in group[name]:
            for r in record:
                if r.name == item:
                    fasta.write(f'>{item}\n{str(r.seq)}\n')
                    break
        fasta.close()
    query_seq.close()

def wether_predict(info_dict):
    '''
    Stop the prediction if there are more than one start/end protein in
    the assemly line
    '''
    start_pro = info_dict['start_pro']
    if len(start_pro) > 1:
        raise Exception(f'There are more than 1 start protein: \
        {",".join(start_pro)}, the protein order cannot be predicted')
    elif len(start_pro) == 1:
        info_dict['start_pro'] = start_pro[0]
    else: 
        info_dict['start_pro'] = ''
        
    end_pro = info_dict['end_pro']
    if len(end_pro) > 1:
        raise Exception(f'There are more than 1 end protein: \
        {",".join(end_pro)}, the protein order cannot be predicted')
    elif len(end_pro) == 1:
        info_dict['end_pro'] = end_pro[0]
    else: 
        info_dict['end_pro'] = ''
        
