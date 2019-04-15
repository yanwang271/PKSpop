#!/usr/bin/env python3

"""
Main script for predicting protein orders in PKS assembly lines.

$ python pks/code/run_analysis.py input_json_file.json

Author: Yan Wang
"""

import os
import subprocess
from sys import argv
import warnings
import json
import extract_seq
import msa
import query_seq
import Ouroboros_analysis
import int_prob
import predict


def parse_input(inpt_path):
    '''
    Parse info from input json file
    '''
    with open(inpt_path) as inpt:
        args = json.load(inpt)

        gbk_path = args['gbk_path']
        pro_id = args['protein_id']
        id_category = args['id_category']
        result_path = args['result_path']
        ouro_path = args['Ouroboros_path']
        ouro_int_frac = args['Ouroboros_int_frac']
        ouro_repeat = args['n_repeat']
    return args


if __name__ == "__main__":

    print('Analysis start......')

    info_dict = parse_input(argv[1])
    result_path = info_dict['result_path']
    if os.path.exists(result_path):
        raise FileExistsError(f"Output path {output_path} already exists")
    else:
        os.mkdir(result_path)
        pred_oupt_path = f'{result_path}/output'
        os.mkdir(pred_oupt_path)

    # Extract whole docking domain sequences from genbank file
    # according to provided protein id information
    extract_seq.extract_seq(info_dict)

    # Cluster the docking domains into 3 classes, then align the
    # class 1 sequences
    msa.clustering(info_dict)
    msa.msa(info_dict)

    # Pair the query sequences and integrate them with the interacting
    # sequences and extra sequences to perform Ouroboros analysis
    query_seq.prepare_query_fl(info_dict)
    print(info_dict)

    # Run Ouroboros analysis with user-defined parameters and find
    # the result with the best LLH
    Ouroboros_analysis.ouroboros_analysis(info_dict)
    print(info_dict)

    # matrix and plot the matrix
    info_dict['output_path'] = pred_oupt_path
    int_prob.prob_mtx(info_dict)
    print(info_dict)

    # Predict the protein order according to interaction probability,
    # start/end protein, protein class
    predict.predict_order(info_dict)
    print(info_dict)
