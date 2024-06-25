from collections import defaultdict
import argparse, os, sys, subprocess
import numpy as np
import time
import pickle as pkl
start = time.time()


#############################################################
########################  Parameters  #######################
#############################################################
parser = argparse.ArgumentParser(description="""PlasGO is a Python library for predicting GO terms of plasmid-encoded proteins.
                                 PlasGO is designed as a hierarchical architecture that leverages the powerful foundation protein language model (ProtTrans-ProtT5) to learn the local context within protein sentences and a BERT model to capture the global context within plasmid sentences.""")
parser.add_argument('--fasta', help='FASTA file of the plasmid-encoded proteins to be annotated, default: example_data/proteins.faa', type=str, default = 'example_data/proteins.faa')
parser.add_argument('--database', help='path of the downloaded pre-annotated database, which consists of 678196 non-redundant plasmid-encoded proteins, along with their GO annotations predicted by PlasGO, default: database', type=str, default='database')
parser.add_argument('--out', help='path to store the annotation results by alignments (Diamond), default: results', type=str, default='results')
parser.add_argument('--threads', help='number of threads utilized for Diamond, default: 2', type=str, default='2')
parser.add_argument('--diamond_mode', help="mode for controlling the sensitivity of Diamond, you can choose one from ['default', 'mid-sensitive', 'sensitive', 'more-sensitive', 'very-sensitive', 'ultra-sensitive'], default: sensitive", type=str, default='sensitive')
inputs = parser.parse_args()


#############################################################
########################  Help info  ########################
#############################################################
def help_info():
    print('')
    print("""Usage of plasgo_diamond.py:
        [--fasta FASTA] FASTA file of the plasmid-encoded proteins to be annotated, default: example_data/proteins.faa
        [--database DATABASE] path of the downloaded pre-annotated database, which consists of 678196 non-redundant plasmid-encoded proteins, along with their GO annotations predicted by PlasGO, default: database
        [--out OUT] path to store the annotation results by alignments (Diamond), default: results
        [--threads THREADS] number of threads utilized for Diamond, default: 2
        [--diamond_mode DIAMOND_MODE] mode for controlling the sensitivity of Diamond, you can choose one from ['default', 'mid-sensitive', 'sensitive', 'more-sensitive', 'very-sensitive', 'ultra-sensitive'], default: sensitive
    """)


#############################################################
######################  Check folders  ######################
#############################################################
db_fn = inputs.database
if not os.path.exists(db_fn):
    print(f"Error! The database folder '{db_fn}' is unavailable. Please use the option '--database' to indicate the directory of the downloaded pre-annotated database folder.")
    help_info()
    sys.exit()


#############################################################
########################  DIAMOND!  #########################
#############################################################
thread_num = inputs.threads
if(inputs.diamond_mode not in ['default', 'mid-sensitive', 'sensitive', 'more-sensitive', 'very-sensitive', 'ultra-sensitive']):
    mode = 'sensitive'
    print(f"Mode selection error! The selected Diamond running mode is not one of ['default', 'mid-sensitive', 'sensitive', 'more-sensitive', 'very-sensitive', 'ultra-sensitive']. Then, running Diamond in sensitive mode ...")
else:
    mode = inputs.diamond_mode
if(mode=='default'):
    mode = ''
else:
    mode = f'--{mode}'
dmnd_path = f'{inputs.database}/database.dmnd'
prot_path = inputs.fasta
res_fn = inputs.out
cmd = f'diamond blastp --threads {thread_num} {mode} -d {dmnd_path} -q {prot_path} -o {res_fn}/resp.tab -k 1 --outfmt 6 qseqid sseqid evalue'
_ = subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


#############################################################
#######################  Write output  ######################
#############################################################
rep_go = pkl.load(open(f'{db_fn}/rep_go.dict', 'rb'))
rep_go = {x:','.join(y) for x,y in rep_go.items()}
f_out = open(f'{res_fn}/results.tsv', 'w')
f_out.write(f'Protein ID\tRepresentative GO term\tGO category\tMember GO terms in the cluster\n')
f=open(f'{res_fn}/resp.tab')
go_dict = pkl.load(open(f'{db_fn}/GO.dict', 'rb'))
ls=f.readlines()
for l in ls:
    d = l.rstrip().split('\t')
    GOs = go_dict[d[1]]
    for Class in ['mf', 'bp', 'cc']:
        if(Class not in GOs):
            continue
        for GO in GOs[Class]:
            f_out.write(f'{d[0]}\t{GO}\t{Class.upper()}\t{rep_go[GO]}\n')
f_out.close()
print(f'The annotation results have been saved in {res_fn}/results.tsv.')


end = time.time()
print(f'The total time for alignment-based annotation against pre-built database is {end-start}s.')
