#!/usr/bin/env python3

"""
Prepare the query sequences file to input into Ouroboros
"""

from Bio import SeqIO
import subprocess

def prepare_query_fl(info_dict):
    '''
    Pair each c-term with all n-term sequence.
    Integrate the paired sequence with the interacting seq pairs 
    and extra seq pairs
    '''
    c_fa = info_dict['c_align_path']
    n_fa = info_dict['n_align_path']
    
    c_pair_fa, n_pair_fa = pairing_seq(c_fa, n_fa)
    
    c_inter_fa = 'PKSpop/data/dd_seq/c_group1.afa'
    c_extra_fa = 'PKSpop/data/dd_seq/c_group1_extra.afa'
    c_query_fl = integrate_fl(c_inter_fa, c_extra_fa, c_pair_fa)
    
    n_inter_fa = 'PKSpop/data/dd_seq/n_group1.afa'
    n_extra_fa = 'PKSpop/data/dd_seq/n_group1_extra.afa'
    n_query_fl = integrate_fl(n_inter_fa, n_extra_fa, n_pair_fa)
    
    info_dict.update({'c_query_path':c_query_fl, 'n_query_path':n_query_fl})
    return info_dict


def pairing_seq(c_fa, n_fa):
    '''
    Pair all c-term sequences to all n_term sequences
    '''
    c_seq = list(SeqIO.parse(c_fa, 'fasta'))
    n_seq = list(SeqIO.parse(n_fa, 'fasta'))
    
    c_pair_fa = c_fa.replace('.afa','_paired.afa')
    c_pair = open(c_pair_fa, 'w')
    
    n_pair_fa = n_fa.replace('.afa','_paired.afa')
    n_pair = open(n_pair_fa, 'w')
    
    for c in c_seq:
        for n in n_seq:
            c_pair.write(f'>{c.name}\n{str(c.seq)}\n')
            n_pair.write(f'>{n.name}\n{str(n.seq)}\n')

    return c_pair_fa, n_pair_fa   
    

def integrate_fl(inter_fa, extra_fa, pair_fa): # all file name
    '''
    Create query sequence file by integrating the query sequence with 
    the interacting seq pairs and extra seq pairs
    '''
    query_fl = f'{pair_fa.replace("paired.afa","ouro_inpt.fasta")}'
    cmd = f'cp {inter_fa} {query_fl}'
    subprocess.check_call(cmd.split(' '))
    
    append_seq(query_fl, extra_fa)
    append_seq(query_fl, pair_fa)
    return query_fl

def append_seq(seq_fl, extra_fl, subset = None):
    '''
    Add the sequence of one file to the other file
    '''
    seq = open(seq_fl,'a')
    if not subset:    
        extra = open(extra_fl).read()
        seq.write(extra)
    else:
        extra = open(extra_fl).readlines()
        for i in subset:
            seq.write(extra[i])
    seq.close()
