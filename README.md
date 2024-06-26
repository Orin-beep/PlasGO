# PlasGO

PlasGO is a Python library __(Linux or Ubuntu only)__ for predicting Gene Ontology (GO) terms of plasmid-encoded proteins. PlasGO is designed as a hierarchical architecture that leverages the powerful foundation protein language model (ProtTrans-ProtT5) to learn the local context within protein sentences and a BERT model to capture the global context within plasmid sentences.

You can use PlasGO in two ways: 1) __alignment-based method__: run the `plasgo_diamond.py` script to align your query proteins against our pre-annotated protein database, which consists of 678,196 non-redundant plasmid-encoded proteins and their predicted high-confidence GO terms; 2) __learning-based method__: run the `plasgo_predict.py` script to predict high-confidence GO terms for your query proteins when they do not have significant alignments to the pre-annotated database.


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


## Download the pre-annotated database (from Zenodo)
Before running the `plasgo_diamond.py` script, you should first download the pre-annotated database of 678,196 plasmid proteins from either [Zenodo](https://zenodo.org/records/12540110/files/database.tar.gz?download=1) or [Google Drive](https://drive.google.com/file/d/1HrWiT_VioxhAoDoK_z6PguHTFDcXQ0gR/view?usp=drive_link).





# Install (Linux or Ubuntu only)
## Dependencies
* [Python 3.x](https://www.python.org/downloads/)
* [NumPy](https://pypi.org/project/numpy/)
* [bidict](https://pypi.org/project/bidict/)
* [Pandas](https://pypi.org/project/pandas/)
* [PyTorch](https://pytorch.org/get-started/previous-versions/)>1.8.0
* [diamond](https://anaconda.org/bioconda/diamond)
* [Prodigal](https://anaconda.org/bioconda/prodigal)
* [biopython](https://pypi.org/project/biopython/)
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)=2.13.0 (conda install -c bioconda blast=2.13.0)
* [HMMER](https://anaconda.org/bioconda/hmmer)
* [treelib](https://pypi.org/project/treelib/)

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


## Format of the output file
The results will be saved in a TSV file (default: `results/host_lineage.tsv`) containing the predicted host lineages from phylum to species level. Each row corresponds to an input plasmid sequence. Examples:

| Contig | Phylum | Class | Order | Family | Genus | Species |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| NZ_CP050042.1  | Pseudomonadota  | Gammaproteobacteria  | Enterobacterales  | Enterobacteriaceae  | Escherichia  | -  | 
| NZ_CP083619.1  | Bacillota  | Clostridia  | Eubacteriales  | Peptostreptococcaceae  | Clostridioides  | Clostridioides difficile  |
| NZ_CP083659.1  | Pseudomonadota  | Gammaproteobacteria  | Moraxellales  | Moraxellaceae  | Acinetobacter  | Acinetobacter variabilis  |
| Z22927.1  | Actinomycetota  | Actinomycetes  | Corynebacteriales  | Corynebacteriaceae  | Corynebacterium  | Corynebacterium glutamicum  |

The dash '-' indicates that HOTSPOT cannot provide accurate predictions for the corresponding input at this taxonomic level.


## Accurate mode with Monte Carlo dropout (MC-dropout)
HOTSPOT provides a specialized *accurate mode* designed to enhance accuracy using MC-dropout based early stop. You can activate this mode using the `--accurate_mode True` option. The results generated in this *accurate mode* will be saved in `results/host_lineage_acc.tsv` by default.


For advanced users, you can utilize the `--mc_num`, `--min_mean`, and `--max_var` options to have control over the specific implementation details of the MC-dropout mechanism. Detailed usage of these options is shown in the following __Full command-line options__ section. Example:

```
python hotspot.py --accurate_mode True
```


## Train your custom models
The data used to train our default models are presented in the table below. If you intend to train custom models, it is necessary to generate corresponding files in the same format, utilizing your own plasmid dataset.

| Required file | Remark | Usage (format) |
| ------------- | ------------- | ------------- |
| [plasmids.fasta](https://drive.google.com/file/d/1LZXSJZCn94JzdzwrCGT1WcMjVT-HWxGx/view?usp=drive_link) | FASTA file of plasmid DNA sequences for training | The downloaded file is a compressed file. You have to unzip it first: `tar -zxvf plasmids.fasta.tar.gz` |
| [host_lineages.tsv](https://drive.google.com/file/d/1a1Tkq_3SywmGo5sVwTC7jzuV0SkfBOIw/view?usp=drive_link) | TSV file containing complete host lineage information from phylum to species for each training plasmid | TSV file with seven columns: plasmid_id, phylum, class, order, family, genus, species |
| [train_val.txt](https://drive.google.com/file/d/15HxvHPkKJx6eXXVNXa88tf-0BJZlBbVp/view?usp=drive_link) | TXT file containing the information of the training/validation sets | The first row should display the list of training plasmids, while the second row should display the list of validation plasmids. Each plasmid in the lists should be separated by a space |

Once you have prepared the files listed in the above table, you can train your custom models by entering into the `train/` folder and sequentially running the `preprocessing_train.py` and `train.py` scripts (in the following case, the three required files are saved in the `HOTSPOT/training_dataset/` folder):
```
cd train/
python preprocessing_train.py --fasta ../training_dataset/plasmids.fasta --host_info ../training_dataset/host_lineages.tsv --train_val_list ../training_dataset/train_val.txt
python train.py
```

## Full command-line options
preprocessing.py:
```
Usage of preprocessing.py:
        [--fasta FASTA] FASTA file of the plasmid DNA sequences to be predicted (either complete sequences or contigs), default: multiple_plasmids.fasta
        [--database DATABASE]   path of the downloaded database folder, which consists of the sequences of PC proteins, MOB/MPF proteins, and replicons, default: database
        [--model_path MODEL_PATH]   path of the folder storing the downloaded or your customized models, default: models
        [--midfolder]   folder to store the intermediate files for prediction, default: temp
        [--len LEN] minimum length of plasmid DNA sequences, default: 1500
        [--threads THREADS] number of threads utilized for preprocessing, default: 2
```

hotspot.py:
```
Usage of hotspot.py:
        [--model_path MODEL_PATH]   path of the folder storing the downloaded or your customized models, default: models
        [--midfolder MIDFOLDER] folder to store the intermediate files (generated by preprocessing.py), default: temp
        [--device DEVICE]   device utilized for prediction ('gpu' or 'cpu'), default: 'gpu'
        [--threads THREADS] number of threads utilized for prediction if 'cpu' is detected ('cuda' not found), default: 2
        [--batch_size BATCH_SIZE]   batch size for prediction, default: 200
        [--out OUT] path to store the prediction results, default: results
        [--accurate_mode ACCURATE_MODE] if this option is set to True, HOTSPOT will run in accurate mode. In accurate mode, the prediction process is slightly slower due to the activation of the MC-dropout mechanism. Additionally, besides the normal output file 'host_lineage.tsv', a supplementary TSV file named 'host_lineage_acc.tsv' will be generated. This file contains high-confidence predictions with a sacrifice of resolution, meaning that for certain inputs, the returned taxa may be at higher levels in the taxonomic hierarchy, default: False
        [--mc_num MC_NUM]   if the accurate mode is activated, you can use this option to specify the number of dropout-enabled forward passes for the MC-dropout mechanism, default: 100
        [--min_mean MIN_MEAN]   the minimum mean value for a prediction that will not trigger early stopping by MC-dropout, default: 0.75
        [--max_var MAX_VAR] the maximum variance value for a prediction that will not trigger early stopping by MC-dropout, default: 0.3
```

preprocessing_train.py:
```
Usage of preprocessing_train.py:
        [--fasta FASTA]    FASTA file of plasmid DNA sequences for training (preferably complete sequences, default: '../training_dataset/plasmids.fasta')
        [--host_info HOST_INFO] TSV file containing complete host lineage information from phylum to species for each training plasmid, default: '../training_dataset/host_lineages.tsv'
        [--database DATABASE]    path of the downloaded database folder, which consists of the sequences of PC proteins, MOB/MPF proteins, and replicons, default: '../database'
        [--model_path MODEL_PATH]  folder to store your customized models, default: 'models'
        [--midfolder MIDFOLDER]   folder to store the intermediate files, default: 'temp'
        [--len LEN] minimum length of plasmid DNA sequences, default: 1500
        [--train_val_list TRAIN_VAL_LIST]  TXT file containing the information of the training/validation sets. The first row should display the list of training plasmids, while the second row should display the list of validation plasmids. Each plasmid in the lists should be separated by a space, default: '../training_dataset/train_val.txt'
        [--train_ratio TRAIN_RATIO] the ratio of the training set size to the total number of input training plasmids. If the train_val_list file is not provided, the training plasmids will be randomly split into train/validation sets using the specified ratio, default: 0.8
        [--labels LABELS]  TXT file containing the specified taxonomic labels for different levels. The file should comprise six rows, where each row corresponds to the labels for the phylum, class, order, family, genus, and species levels, respectively. Within each row, the labels should be separated by tabs ('	'), default: None
        [--num_plasmids NUM_PLASMIDS]    minimun number of training plasmids associated with a taxonomic label (used when the labels file is not provided), default: 20
        [--add_frags ADD_FRAGS]   whether to augment the training set by randomly cutting fragments ranging from 1.5 to 15 kbp (this may slightly slow down the training process, but it will significantly enhance the performance of host prediction), default: False
        [--num_frags NUM_FRAGS]   maximum number of added fragments from each training plasmid, default: 5
        [--threads THREADS] number of threads utilized for preprocessing, default: 2
```

train.py
```
Usage of train.py:
        [--model_path MODEL_PATH]   folder to store your customized models, default: models
        [--midfolder MIDFOLDER] folder to store the intermediate files (generated by preprocessing_train.py), default: temp
        [--device DEVICE]   device utilized for training ('gpu' or 'cpu'), default: 'gpu'
        [--threads THREADS] number of threads utilized for training if 'cpu' is detected ('cuda' not found), default: 2
        [--lr LR]   learning rate for training the models, default: 0.005
        [--batch_size BATCH_SIZE]   batch size for training the models, default: 200
        [--epoch_num EPOCH_NUM] number of epochs for training the models, default: 20
        [--dropout DROPOUT] dropout rate for training the models, default: 0.5 
```


## References
How to cite this tool:
```
Ji, Yongxin, et al. "HOTSPOT: hierarchical host prediction for assembled plasmid contigs with transformer." Bioinformatics 39.5 (2023): btad283. https://academic.oup.com/bioinformatics/article/39/5/btad283/7136643
```
