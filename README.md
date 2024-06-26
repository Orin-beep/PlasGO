# PlasGO

PlasGO is a Python library __(Linux or Ubuntu only)__ for predicting Gene Ontology (GO) terms of plasmid-encoded proteins. PlasGO is designed as a hierarchical architecture that leverages the powerful foundation protein language model (ProtTrans-ProtT5) to learn the local context within protein sentences and a BERT model to capture the global context within plasmid sentences.

You can use PlasGO in two ways: 1) __alignment-based method__: run the `plasgo_diamond.py` script to align your query proteins against our pre-annotated protein database, which consists of 678,196 non-redundant plasmid-encoded proteins and their predicted high-confidence GO terms by PlasGO; 2) __learning-based method__: run the `plasgo_predict.py` script to predict high-confidence GO terms for your query proteins when they do not have significant alignments to the pre-annotated database.


### E-mail: yongxinji2-c@my.cityu.edu.hk


### Version: V1.0 (2024-06-26)


# 1. Alignment-based Method (using Diamond)
## Install
If you only need to annotate your plasmid proteins by aligning them against the pre-annotated database using Diamond, all you require is an environment with [Python 3.x](https://www.python.org/downloads/).
```
git clone https://github.com/Orin-beep/PlasGO
tar -jxvf diamond.tar.bz2
chmod 755 diamond
```


## Download the pre-annotated database
Before running the `plasgo_diamond.py` script, you should first manually download the pre-annotated database of 678,196 plasmid proteins from either [Zenodo](https://zenodo.org/records/12540110/files/database.tar.gz?download=1) or [Google Drive](https://drive.google.com/file/d/1HrWiT_VioxhAoDoK_z6PguHTFDcXQ0gR/view?usp=drive_link).

After downloading the file `models.tar.gz`, place it in the same directory as `plasgo_diamond.py` and uncompress it:
```
tar -zxvf database.tar.gz
rm database.tar.gz
```


## Usage
Then, you can easily run `plasgo_diamond.py` to annotate your query plasmid proteins using Diamond. Specifically, the query proteins will be annotated with the same GO terms as the target proteins with the best hits in the database.


### Simple example
```
python plasgo_diamond.py --fasta example_data/proteins.faa --database database --threads 8
```
The annotation results will be saved in `results/results.tsv` by default.


# 2. Learning-based Method (predict using the PlasGO models)
## Install (Linux or Ubuntu only)
### Dependencies
* [Python 3.x](https://www.python.org/downloads/)
* [NumPy](https://pypi.org/project/numpy/)
* [PyTorch](https://pytorch.org/get-started/previous-versions/)>1.8.0
* [biopython](https://pypi.org/project/biopython/) (pip install biopython)
* [datasets](https://pypi.org/project/datasets/) (pip install datasets)
* [transformers](https://huggingface.co/docs/transformers/installation) (pip install transformers[sentencepiece])
* [sentencepiece](https://pypi.org/project/sentencepiece/) (pip install sentencepiece)

If you want to use GPU to accelerate the program:
* CUDA
* PyTorch-GPU
* For CPU version PyTorch: ```conda install pytorch torchvision torchaudio cpuonly -c pytorch```
* For GPU version PyTorch: search [PyTorch](https://pytorch.org/get-started/previous-versions/) to find the correct CUDA version according to your computer


## Prepare the environment
After cloning this repository (```git clone https://github.com/Orin-beep/HOTSPOT```), you can use Anaconda to install ```environment.yaml```. This will install all packages you need in GPU mode (make sure you have installed CUDA on your system to use the GPU version; otherwise, HOTSPOT will run in CPU mode). The installing command is: 
```
git clone https://github.com/Orin-beep/HOTSPOT
cd HOTSPOT/
conda env create -f environment.yaml -n hotspot
conda activate hotspot
```
If Anaconda fails to work, you can prepare the environment by individually installing the packages listed in the __Dependencies__ section.

## Prepare default database and models (from Google Drive)
To download the default database and models, you can use the following bash scripts: 
```
sh prepare_db.sh        # download and unzip the database folder, 91.8 MB
sh prepare_mdl.sh       # download and unzip the model folder, 1.78 GB
```


## Alternative way: download default database and models manually
If the above bash scripts do not work, you can manually download the default database and models using the following Google Drive links:
* [database.tar.gz](https://drive.google.com/file/d/1ZSTz3kotwF8Zugz_aBGDtmly8BVo9G4T/view)
* [models.tar.gz](https://drive.google.com/file/d/1bnA1osvYDgYBi-DRFkP-HrvcnvBvbipF/view)

After downloading the `database.tar.gz` and `models.tar.gz` to HOTSPOT's main directory, you have to unzip them:
```
tar -zxvf database.tar.gz
rm database.tar.gz

tar -zxvf models.tar.gz
rm models.tar.gz
```


# Usage
Before running HOTSPOT, you should run `preprocessing.py` to encode the input plasmid sequences into sentences. After that, you can use `hotspot.py` for host prediction.

## Simple example
```
python preprocessing.py --fasta example_plasmids/NZ_CP083659.fasta --database database --model_path models
python hotspot.py
```


# Format of the output file
The annotation results will be saved in a TSV file (default: `results/results.tsv`) containing the multiple cluster-level GO term labels annotated to each query plasmid protein. An output example of a query protein `WP_071931462`:

| Protein ID | Representative GO term | GO category | Probability (only available for `plasgo_predict.py`) | Member GO terms in the cluster |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| WP_071931462  | GO:0003723  | MF  | 0.5  | GO:0003677,GO:0003723  |
| WP_071931462  | GO:1901363  | MF  | 1.0  | GO:1901363  |
| WP_071931462  | GO:0003676  | MF  | 1.0  | GO:0003676  |
| WP_071931462  | GO:0097159  | MF  | 1.0  | GO:0097159  |
| WP_071931462  | GO:0005488  | MF  | 1.0  | GO:0005488  |
| WP_071931462  | GO:0009987  | BP  | 0.5  | GO:0009987  |
| WP_071931462  | GO:0016020  | CC  | 0.99  | GO:0016020  |
| WP_071931462  | GO:0110165  | CC  | 1.0  | GO:0110165  |

The second column shows the representative GO terms of the predicted clusters, as determined by the REVIGO tool. The fourth column lists the member GO terms within each cluster. Importantly, if a query protein is predicted to be associated with a particular cluster, it means the protein is annotated to at least one of the GO terms within that cluster.


# Full command-line options
plasgo_diamond.py:
```
Usage of plasgo_diamond.py:
        [--fasta FASTA] FASTA file of the plasmid-encoded proteins to be annotated, default: example_data/proteins.faa
        [--database DATABASE] path of the downloaded pre-annotated database, which consists of 678196 non-redundant plasmid-encoded proteins, along with their GO annotations predicted by PlasGO, default: database
        [--out OUT] path to store the annotation results by alignments (Diamond), default: results
        [--threads THREADS] number of threads utilized for Diamond, default: 2
        [--diamond_mode DIAMOND_MODE] mode for controlling the sensitivity of Diamond, you can choose one from ['default', 'mid-sensitive', 'sensitive', 'more-sensitive', 'very-sensitive', 'ultra-sensitive'], default: sensitive
```

preprocessing.py:
```
Usage of preprocessing.py:
        [--fasta FASTA] FASTA file of the plasmid-encoded proteins to be predicted, default: example_data/proteins.faa
        [--plasmids PLASMIDS] TXT file containing the contextual information if available (not mandatory), where the proteins are arranged in the same order as their encoding within the plasmid. Each row includes the protein IDs separated by semicolon (';') within one complete plasmid or plasmid segment. For the proteins without contextual information, PlasGO will predict them mainly using the embedding layers of the models. The example data offers three complete plasmids represented by the IDs of their encoded proteins. default: example_data/3plasmids.txt
        [--midfolder MIDFOLDER] folder to store the intermediate files for prediction, default: temp
        [--prott5 PROTT5] folder to store the downloaded ProtT5-XL-U50 model, default: prot_t5_xl_uniref50
        [--device DEVICE] device utilized for generating original per-protein embeddings with ProtT5 ('gpu' or 'cpu'), default: 'gpu'
        [--batch_size BATCH_SIZE] batch size for protein embedding with the ProtT5 model. If your GPU is out of memory, you can try to reduce this parameter to 1, default: 8
```

plasgo_predict.py:
```
Usage of plasgo_predict.py:
        [--midfolder MIDFOLDER] folder to store the intermediate files (generated by preprocessing.py), default: temp
        [--model_path MODEL_PATH] path of the folder storing the downloaded or your customized models, default: models
        [--out OUT] path to store the prediction results, default: results
        [--min_prob MIN_PROB] the minimum probability for determining a high-confidence GO term prediction (the maximum value is limited to 0.3), default: 0.3
        [--device DEVICE] device utilized for GO term prediction ('gpu' or 'cpu'), default: 'gpu'
        [--batch_size BATCH_SIZE] batch size (plasmid sentence count in a batch) for GO term prediction. If your GPU is out of memory during prediction, you can try to reduce this parameter, default: 32
```
