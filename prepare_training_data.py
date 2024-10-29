import pickle as pkl
from Bio import SeqIO
import os, argparse
import csv
from collections import defaultdict
import numpy as np


# This script is used to generate dataset in the format which can be directly input to the PlasGO program (either training or prediction)
#############################################################
########################  Parameters  #######################
#############################################################
parser = argparse.ArgumentParser(description="""PlasGO is a Python library for predicting GO terms of plasmid-encoded proteins.
                                 PlasGO is designed as a hierarchical architecture that leverages the powerful foundation protein language model (ProtTrans-ProtT5) to learn the local context within protein sentences and a BERT model to capture the global context within plasmid sentences.""") 
parser.add_argument('--input_folder', help='folder storing raw data for processing the dataset, default: raw_data', type=str, default = 'raw_data')
parser.add_argument('--out', help='folder to store the processed dataset for training, default: training_data', type=str, default='training_data')
inputs = parser.parse_args()


#############################################################
########################  Help info  ########################
#############################################################
def help_info():
    print('')
    print("""Usage of prepare_training_data.py:
        [--input_folder INPUT_FOLDER]   folder storing raw data for processing the dataset, default: raw_data
        [--out OUT] folder to store the processed dataset for training, default: training_data
    """)


out_fn = inputs.out
if not os.path.isdir(out_fn):
    os.makedirs(out_fn)


# Plasmids containing more than 56 proteins were divided into multiple segments with an overlap of 14 proteins (1/4 of the maximum length).
def segmentation(pls_prot):
    length = 56 # hardcode, length 56, the median of all plasmids' lengths
    stride = int(length*0.75)

    corpus = []
    for pls in pls_prot:
        sentence = pls_prot[pls]

        init = 0
        while(init+length<=len(sentence)):
            segment = sentence[init:init+length]
            corpus.append(segment)
            init+=stride
        # the last segment
        if(len(sentence)<length or (len(sentence)-length)%stride!=0):
            segment = sentence[-length:]
            corpus.append(segment)
    print(f'Number of cut segments before deduplication:', len(corpus))

    # deduplication
    new_corpus = set()
    for l in corpus:
        string = '|'.join(l)
        new_corpus.add(string)
    print('Number of segments after deduplication:', len(new_corpus))

    # decode
    corpus = []
    for string in new_corpus:
        d = string.split('|')
        corpus.append(d)
    return corpus


# write into csv files which can be read by the datasets module
def prepare_datasets(out_fn, corpus, label_dict):
    length = 56
    header = ['labels', 'proteins']
    f = open(out_fn, 'w')
    writer = csv.writer(f)
    writer.writerow(header)
    none_label = 'None'
    prots = set(label_dict)

    for List in corpus:
        if(prots&set(List)==set()):
            continue

        row = []
        labels, proteins = [], []

        for prot in List:
            # get protein IDs and labels
            if(prot not in prots):  # do not have labels
                labels.append(none_label)
                proteins.append(f'*{prot}')
            else:
                label = label_dict[prot]
                label = [str(i) for i in label]
                label = ' '.join(label)
                labels.append(label)
                proteins.append(prot)

        while(len(labels)!=length):
            labels.append(none_label)
            proteins.append(none_label)

        labels = ';'.join(labels)
        proteins = ';'.join(proteins)
        row.append(labels)
        row.append(proteins)
        writer.writerow(row)
    f.close()


def prepare_plasgo_dataset(pls_prot, fn, data_splitting, label_encodings):
    corpus = segmentation(pls_prot)
    
    for Class in ['mf', 'bp', 'cc']:    # three categories
        out_fn = f'{fn}/dataset_{Class}'
        if not os.path.exists(out_fn):
            os.mkdir(out_fn)
        true_label = label_encodings[Class]
        train, val, test = data_splitting[Class][0], data_splitting[Class][1], data_splitting[Class][2]
        
        train_dict = {x:y for x,y in true_label.items() if x in train}
        val_dict = {x:y for x,y in true_label.items() if x in val}
        test_dict = {x:y for x,y in true_label.items() if x in test}
        
        prepare_datasets(f'{out_fn}/training.csv', corpus, train_dict)
        prepare_datasets(f'{out_fn}/validation.csv', corpus, val_dict)
        prepare_datasets(f'{out_fn}/test.csv', corpus, test_dict)


def check_consistent(path, variable):
    X = pkl.load(open(path, 'rb'))
    print(X==variable)


def faa2set(path):
    res = set()
    for s in SeqIO.parse(path, 'fasta'):
        res.add(s.id)
    return res


def load_input(fn, out_fn):
    class_name = {'cc': 'cellular_component', 'mf': 'molecular_function', 'bp': 'biological_process'}
    
    # labels
    label_file = f'{fn}/Label_info.tsv'
    labels = defaultdict(dict)
    rep2go = {}
    with open(label_file, mode='r', newline='', encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            labels[row['Category'].lower()][row['Label']] = int(row['Label_index'])
            rep2go[row['Label']] = set(row['Member_GO_terms'].split(','))
    pkl.dump(labels, open(f'{out_fn}/labels.dict', 'wb'))
    pkl.dump(rep2go, open(f'{out_fn}/rep_go.dict', 'wb'))

    # pls_prot
    plasmid_file = f'{fn}/Plasmid_info.tsv'
    pls_prot = defaultdict(list)
    with open(plasmid_file, mode='r', newline='', encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            pls_prot[row['Plasmid_ID']] = row['Order_Encoded_proteins'].split(',') 

    # label encoding
    label_encodings = {}
    for Class in ['bp', 'cc', 'mf']:
        prot_label_path = f'{fn}/{class_name[Class]}/Protein_labels.tsv'
        prot_labels = {}
        with open(prot_label_path, mode='r', newline='', encoding='utf-8') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                prot_labels[row['Protein_ID']] = set(row['Annotated_labels'].split(','))
        labels_tmp = labels[Class]
        label_encodings_tmp = {}
        for prot in prot_labels:
            tmp = np.zeros(len(labels_tmp))
            for label in prot_labels[prot]:
                idx = labels_tmp[label]
                tmp[idx] = 1
            label_encodings_tmp[prot] = tmp
        label_encodings[Class] = label_encodings_tmp

    # data_splitting
    data_splitting = {}
    for Class in ['bp', 'cc', 'mf']:
        path = f'{fn}/{class_name[Class]}/'
        train = faa2set(f'{path}/training.faa')
        validation = faa2set(f'{path}/validation.faa')
        test = faa2set(f'{path}/test.faa')
        data_splitting_tmp = [train, validation, test]
        data_splitting[Class] = data_splitting_tmp
    return pls_prot, data_splitting, label_encodings


input_folder = inputs.input_folder
pls_prot, data_splitting, label_encodings = load_input(input_folder, out_fn)
prepare_plasgo_dataset(pls_prot, out_fn, data_splitting, label_encodings)


