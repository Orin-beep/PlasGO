import argparse, os, sys
import re
import numpy as np
import time
import torch
import pickle as pkl
from transformers import T5Tokenizer, T5EncoderModel
from Bio import SeqIO
start = time.time()


#############################################################
########################  Parameters  #######################
#############################################################
parser = argparse.ArgumentParser(description="""PlasGO is a Python library for predicting GO terms of plasmid-encoded proteins.
                                 PlasGO is designed as a hierarchical architecture that leverages the powerful foundation protein language model (ProtTrans-ProtT5) to learn the local context within protein sentences and a BERT model to capture the global context within plasmid sentences.""") 
parser.add_argument('--fasta', help='FASTA file of the plasmid-encoded proteins to be predicted, default: example_data/proteins.faa', type=str, default = 'example_data/proteins.faa')
parser.add_argument('--plasmids', help="TXT file containing the contextual information if available (not mandatory), where the proteins are arranged in the same order as their encoding within the plasmid. Each row includes the protein IDs separated by semicolon (';') within one complete plasmid or plasmid segment. For the proteins without contextual information, PlasGO will predict them mainly using the embedding layers of the models. The example data offers three complete plasmids represented by the IDs of their encoded proteins. default: example_data/3plasmids.txt", type=str, default = 'example_data/3plasmids.txt')
parser.add_argument('--midfolder', help='folder to store the intermediate files for prediction, default: temp', type=str, default='temp')
parser.add_argument('--prott5', help='folder to store the downloaded ProtT5-XL-U50 model, default: prot_t5_xl_uniref50', type=str, default='prot_t5_xl_uniref50')
parser.add_argument('--device', help="device utilized for generating original per-protein embeddings with ProtT5 ('gpu' or 'cpu'), default: 'gpu'", type=str, default = 'gpu')
parser.add_argument('--batch_size', help="batch size for protein embedding with the ProtT5 model. If your GPU is out of memory, you can try to reduce this parameter to 1, default: 8", type=int, default = 8)
inputs = parser.parse_args()


#############################################################
########################  Help info  ########################
#############################################################
def help_info():
    print('')
    print("""Usage of preprocessing.py:
        [--fasta FASTA] FASTA file of the plasmid-encoded proteins to be predicted, default: example_data/proteins.faa
        [--plasmids PLASMIDS] TXT file containing the contextual information if available (not mandatory), where the proteins are arranged in the same order as their encoding within the plasmid. Each row includes the protein IDs separated by semicolon (';') within one complete plasmid or plasmid segment. For the proteins without contextual information, PlasGO will predict them mainly using the embedding layers of the models. The example data offers three complete plasmids represented by the IDs of their encoded proteins. default: example_data/3plasmids.txt
        [--midfolder MIDFOLDER] folder to store the intermediate files for prediction, default: temp
        [--prott5 PROTT5] folder to store the downloaded ProtT5-XL-U50 model, default: prot_t5_xl_uniref50
        [--device DEVICE] device utilized for generating original per-protein embeddings with ProtT5 ('gpu' or 'cpu'), default: 'gpu'
        [--batch_size BATCH_SIZE] batch size for protein embedding with the ProtT5 model. If your GPU is out of memory, you can try to reduce this parameter to 1, default: 8
    """)


#############################################################
######################  Check folders  ######################
#############################################################
out_fn = inputs.midfolder
if not os.path.isdir(out_fn):
    os.makedirs(out_fn)

prot_path = inputs.fasta
plasmid_path = inputs.plasmids
prott5_folder = inputs.prott5
if not os.path.exists(prott5_folder):
    print(f"Error! The ProtT5-XL-U50 folder '{prott5_folder}' is unavailable. Please download the ProtT5-XL-U50 model using the command 'wget https://zenodo.org/record/4644188/files/prot_t5_xl_uniref50.zip' and unzip it. Then, you can use the option '--prott5' to indicate the directory of the downloaded ProtT5-XL-U50 model.")
    help_info()
    sys.exit()


#############################################################
####################  Write dataset file  ###################
#############################################################
csv_fn = f'{out_fn}/plasmids.csv'

# read FASTA file
prot_idx = {}
prot_id_set = set()
idx = 1
all_sequences = []
for s in SeqIO.parse(prot_path, 'fasta'):
    if(s.id in prot_id_set):
        print(f'Repeated protein {s.id} detected! PlasGO will use the sequence of the first {s.id} only.')
        continue
    prot_id_set.add(s.id)
    all_sequences.append(str(s.seq))
    prot_idx[s.id] = idx
    idx+=1
pkl.dump(prot_idx, open(f'{out_fn}/prot_idx.dict', 'wb'))

# generate corpus
corpus = []
length = 56
stride = int(length*0.75)
if not os.path.exists(plasmid_path):
    print(f'Contextual information is not provided. PlasGO will predict GO terms of the input proteins individually.')
    for prot in prot_id_set:
        corpus.append([prot])
else:
    f=open(plasmid_path)
    ls=f.readlines()
    for l in ls:
        raw_sentence = l.rstrip().split(';')
        sentence = []
        for prot in raw_sentence:
            if(prot not in prot_id_set):
                print(f'Warning! Protein {prot} is not included in the input FASTA file {prot_path} and will not be predicted.')
            else:
                sentence.append(prot)
        if(len(sentence)<=length):
            corpus.append(sentence)
        else:
            init = 0
            while(init+length<=len(sentence)):
                segment = sentence[init:init+length]
                corpus.append(segment)
                init+=stride
            if((len(sentence)-length)%stride!=0):
                segment = sentence[-length:]
                corpus.append(segment)
    # deduplication
    new_corpus = set()
    for l in corpus:
        string = '|'.join(l)
        new_corpus.add(string)
    # decode
    corpus = []
    for string in new_corpus:
        d = string.split('|')
        corpus.append(d)

    # check singletons
    inuse_prots = set()
    for sentence in corpus:
        for prot in sentence:
            inuse_prots.add(prot)
    singletons = prot_id_set - inuse_prots
    for prot in singletons:
        corpus.append([prot])


# write CSV
f=open(csv_fn, 'w')
f.write(f'proteins\n')
for sentence in corpus:
    sentence = sentence+['None']*(length-len(sentence))
    row = ';'.join(sentence)
    f.write(f'{row}\n')
f.close()


#############################################################
###############  Learn embeddings with ProtT5  ##############
#############################################################
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
device_opt = inputs.device
if(device==torch.device("cuda")):
    print("GPU detected. Running protein embedding with GPU ...")
else:
    if(device_opt=='cpu'):
        print("GPU not detected. Running protein embedding with CPU ...")
    else:
        print("GPU not detected. We highly recommend you to run this script with GPU because of the large size of the ProtT5 model. If you still want to run with CPU, please specify the option '--device cpu'")

tokenizer = T5Tokenizer.from_pretrained(prott5_folder, do_lower_case=False)
model = T5EncoderModel.from_pretrained(prott5_folder).to(device)
batch_size = inputs.batch_size
batch_num = int(len(all_sequences)/batch_size)+1
features = []
for idx in range(batch_num):
    sequences = all_sequences[batch_size*idx:batch_size*(idx+1)]
    sequences = [" ".join(list(re.sub(r"[UZOB]", "X", sequence))) for sequence in sequences]
    if(sequences==[]):
        continue
    ids = tokenizer.batch_encode_plus(sequences, add_special_tokens=True, padding=True)
    input_ids = torch.tensor(ids['input_ids']).to(device)
    attention_mask = torch.tensor(ids['attention_mask']).to(device)
    with torch.no_grad():
        embedding = model(input_ids=input_ids,attention_mask=attention_mask)
    embedding = embedding.last_hidden_state
    for seq_num in range(len(sequences)):   # Global Average Pooling (GAP)
        seq_len = (attention_mask[seq_num] == 1).sum()
        seq_emd = embedding[seq_num][:seq_len - 1]
        emd_per_protein = seq_emd.mean(dim=0)
        features.append(emd_per_protein.cpu().numpy())

features = np.vstack(features)
print(f'Protein embedding finished! The shape of the embedded matrix is {features.shape}.')
pkl.dump(features, open(f'{out_fn}/raw_embeddings.pkl', 'wb'))


end = time.time()
print(f'The total time for preprocessing (corpus generation + protein embedding) is {end-start}s.')
