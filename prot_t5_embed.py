from transformers import T5Tokenizer, T5EncoderModel
from Bio import SeqIO
import argparse, os
import pickle as pkl
import numpy as np
import re
import torch


# The code for extracting per-protein embeddings using PLM 'prot_t5_xl_uniref50'.
# For more details of the model 'prot_t5_xl_uniref50', please refer to https://huggingface.co/Rostlab/prot_t5_xl_uniref50.
# You only need to input the protein sequences to be embedded in FASTA format.


#############################################################
########################  Parameters  #######################
#############################################################
parser = argparse.ArgumentParser(description="The script for extracting per-protein embeddings using PLM 'prot_t5_xl_uniref50'.")
parser.add_argument('--faa', help='path of the protein FASTA file to be embedded, default: raw_data/proteins.faa', type=str, default='raw_data/proteins.faa')
parser.add_argument('--out', help='folder to store the protein embedding results, default: training_data/protein_embeddings', type=str, default='training_data/protein_embeddings/')
parser.add_argument('--device', help="device utilized for protein embedding ('gpu' or 'cpu'), default: 'gpu'", type=str, default = 'gpu')
parser.add_argument('--batch_size', help="batch size used for protein embedding. If your GPU is out of memory, you can try to reduce this parameter, default: 16", type=int, default = 16)
parser.add_argument('--prott5', help='folder to store the downloaded ProtT5-XL-U50 model, default: prot_t5_xl_uniref50/', type=str, default='prot_t5_xl_uniref50/')
inputs = parser.parse_args()


#############################################################
########################  Help info  ########################
#############################################################
def help_info():
    print('')
    print("""Usage of prot_t5_embed.py:
        [--faa FAA] path of the protein FASTA file to be embedded, default: raw_data/proteins.faa
        [--out OUT] folder to store the protein embedding results, default: training_data/protein_embeddings
        [--device DEVICE]   device utilized for protein embedding ('gpu' or 'cpu'), default: 'gpu'
        [--batch_size BATCH_SIZE]   batch size used for protein embedding. If your GPU is out of memory, you can try to reduce this parameter, default: 16
        [--prott5 PROTT5]   folder to store the downloaded ProtT5-XL-U50 model, default: prot_t5_xl_uniref50/
    """)


out_fn = inputs.out
if not os.path.exists(out_fn):
    os.mkdir(out_fn)


# read protein sequences
seq_dict = {}
for s in SeqIO.parse(inputs.faa, 'fasta'): 
    seq_dict[s.id] = str(s.seq)

prot_list = sorted(list(seq_dict))
prot2idx = {}
all_sequences = []
idx = 1
for prot in prot_list:
    prot2idx[prot] = idx
    all_sequences.append(seq_dict[prot])
    idx+=1
pkl.dump(prot2idx, open(f'{out_fn}/prot2idx.dict', 'wb'))


# device
if(inputs.device=='cpu'):
    device = torch.device("cpu")
    print("Running protein embedding with CPU ...")
else:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if(device==torch.device("cuda")):
        print("GPU detected. Running protein embedding with GPU ...")        
    else:
        print("GPU not detected. Running protein embedding with CPU ...")


# load model
tokenizer = T5Tokenizer.from_pretrained(inputs.prott5, do_lower_case=False)
model = T5EncoderModel.from_pretrained(inputs.prott5).to(device)


# run!
batch_size = inputs.batch_size
print(f'Run protein embedding with batch size of {batch_size} ...')
batch_num = int(len(all_sequences)/batch_size)+1
features = []
for idx in range(batch_num):
    sequences = all_sequences[batch_size*idx:batch_size*(idx+1)]
    sequences = [" ".join(list(re.sub(r"[UZOB]", "X", sequence))) for sequence in sequences]
    ids = tokenizer.batch_encode_plus(sequences, add_special_tokens=True, padding=True)
    input_ids = torch.tensor(ids['input_ids']).to(device)
    attention_mask = torch.tensor(ids['attention_mask']).to(device)
    
    with torch.no_grad():
        embedding = model(input_ids=input_ids,attention_mask=attention_mask)
    
    embedding = embedding.last_hidden_state.cpu()
    for seq_num in range(len(sequences)):
        seq_len = (attention_mask[seq_num] == 1).sum()
        seq_emd = embedding[seq_num][:seq_len - 1]
        emd_per_protein = seq_emd.mean(dim=0)
        features.append(emd_per_protein.numpy())


# save into file
features = np.vstack(features)
pkl.dump(features, open(f'{out_fn}/prot_embeddings.pkl', 'wb'))
print(f'Protein embedding finished. The feature matrix of shape {features.shape} have been saved as {out_fn}/prot_embeddings.pkl.')
