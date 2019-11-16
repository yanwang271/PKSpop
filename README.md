PKSpop
==============

A pipeline to predict type I PKSs protein order in polyketide biosynthetic assembly lines.


PKSpop comprises three main steps to infer protein order: 
1. Identify class memberships for query docking domains and align the sequences
2. Pair each class I Ndd with all class I Cdds and use Ouroboros to predict the interaction probability for each pair. The probabilities are filled into a matrix
3. Infer protein order from the probability matrix by the Hungarian algorithm, a global optimization method, which finds a match between Ndds and Cdds that maximizes the overall interaction probability. The inference method takes the assembly line constraints and compatibility class into account.

## Run PKSpop
```
python PKSpop/code/run_analysis.py input_json_file.json
```

## Input
The input of PKSpop is a JSON file which contains following informaion:
* gbk_path: path where the antiSMASH ```.gbk``` file of the query PKS gene cluster
* protein_id: list of identifiers of query proteins whose order will be predicted
* id_category: the category of the protein identifiers: "gene", "protein_id" or "locus_tag"
* result_path: path where result will be saved
* Ouroboros_path: path to Ouroboros repositry
* Ouroboros_int_frac: list of ```int_frac``` parameter of Ouroboros with default ```[0.9, 0.8]```. It is recommended to add numbers below 0.8 if there are more than 10 query proteins.
* n_repeat: number of repeat time to run Ouroboros with each ```int_frac``` parameter
An example can be found in ```data/test```

## Output
The prediction results are in ```result_path/output/```
* protein_order_prediction.txt gives the predicted order of the query assembly line
* protein_interaction_prediction.txt gives the predicted interacting protein pairs in the query assembly line
* int_prob_mtx.csv is the pairwise interaction probabilities matrix of query proteins predicted by Ouroboros

Additional files used in prediction process:
* dd_raw.fasta is the raw sequences extracted from the input ```.gbk``` file
* dd_class_*.fasta is the sequences of 3 compatibility classes
* dd_hmmscan_oupt.txt is the result of hmmscan, which contains the class information of the sequences
* dd_class_1_aln.afa is the aligned class 1 sequences
* dd_class_1.afa is the conserved region on the aligned sequences
* dd_class_1_paired.afa is the paired sequences of all query proteins
* dd_class_1_ouro_inpt.fasta is the fasta file that input into Ouroboros
* Ouroboros_class_1_ouro_inpt_*_soft_warm_* is the Ouroboros' output 

## License
This project is licensed under the BSD-3 license. See the LICENSE file for details.

## Requirements
PKSpop requires Python 3.6+. The following tools should be installed/downloaded before running PKSpop:
* [Ouroboros](https://github.com/miguelcorrea/Ouroboros)
* [HMMER](https://hmmer.org)
Packages:
* Biopython
* NumPy
* Pandas
* SciPy


