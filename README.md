# pep2prot
pep2prot (DIA-NN edition) maps peptide sequences to proteins when annotation was not provided by the used proteomics software (here DIA-NN) and returns info about the tryptic state of each peptide.

## Requirements

- Input: path to a DIA-NN output csv file, semi-colon separated, with unique peptide sequences in a "Peptide" column

- Input: path to the proteome database in FASTA format which was used for the DIA-NN search

- path to the output folder

## Setup

Create a conda environment (use [Miniconda](https://docs.conda.io/en/latest/miniconda.html)):

`conda create -n pep2prot`

`conda activate pep2prot`

Install NumPy, Pandas, PyYAML, BioPython and termcolor:

`conda install -c conda-forge -c bioconda numpy pandas pyyaml biopython termcolor`

Clone this repo:

`git clone https://github.com/npinter/pep2prot.git`

## How To Use

1. Open the `/config` folder, copy and paste the YAML file ending with `.diann.yaml`.
2. Edit the copied YAML file with any editor and change paths.
```
num_jobs: number of parallel jobs, 0=all availiable threads
protein_expression_csv: output csv from DIA-NN
fasta_db_path: FASTA database used in the DIA-NN search
output_folder: output folder for pep2prot result
```
3. Save the YAML file and start `pep2prot_diann.bat` (Windows, modify the batch file first) or `pep2prot_diann.py`.
4. Select your YAML file in the prompt.
5. Press `Enter` to start the mapping.
