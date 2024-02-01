Subsampling Tool

Authors: Gerlinde Grabmann, Christina Pfeiffer

# Overview

This Python script is designed to parse (compressed) FASTA or FASTQ files and randomly select sequences to generate a subsample.
It provides options for fast parsing, preserving the original formatting of sequences, and setting a seed for reproducibility.

# Usage

## Requirements
    Python 3.9.18
    Biopython library

## Repository
    The script, YAML file, and example files can be found in the following repo:
    <https://github.com/tinscherl2/Subsample_biological_sequences>

## Install the required dependencies:
    <conda install -c conda-forge biopython> (see https://biopython.org/wiki/Packages)
    Alternatively you can create a conda environment with the required dependencies from the environment.yml file provided:
    <conda create --name <env> --file environment.yml>

## Running the Script
Use the script from the command line with the following format:

python final_Grabmann_Pfeiffer.py INPUT_FILE NUM_READS [-f] [-p] [-s SEED]
or
python final_Grabmann_Pfeiffer.py INPUT_FILE NUM_READS [--fast] [--preserve-format] [--seed SEED]

INPUT_FILE: Path to the input file (FASTA or FASTQ, can be gzipped).
NUM_READS: Number of reads to subsample.
-f or --fast: Use fast parsing. Do not use option with files that do not fit into memory. 
-p or --preserve-format: Preserve the original formatting of sequences. The number of lines will not change. This option is not compatible with with the fast flag.
-s SEED or --seed SEED: Seed for the random number generator for reproducibility or paired-end reads. Use the same seed for the forward and the reverse read file. Replace SEED by your favorite number.

## Examples 

- Subsample 3 reads from a FASTA file
python final_Grabmann_Pfeiffer.py.py fasta_testfile.fasta 3

- Subsample 3 reads from a FASTA file with fast parsing and setting seed
python final_Grabmann_Pfeiffer.py.py fasta_testfile.fasta 3 -f -s 42

- Subsample 3 reads from a FASTA file preserving the format and setting seed 42
python final_Grabmann_Pfeiffer.py.py fasta_testfile.fasta 3 -p -s 42


# Notes

The script supports both FASTA and FASTQ file formats, and it can handle gzipped files.

Importantly, large files that do not fit into memory should not use the fast option. 
Original format cannot be preserved for fast parsing.

The -f flag enables fast parsing for large files that fit into memory.
The -p flag preserves the original formatting of sequences.
The -s flag allows setting a seed for reproducibility.

The fasta_festfile.fasta was generated with ChatGPT.
The paired-end fastq.gz files were obtained from https://zenodo.org/records/1236641.
The faulty.fasta is a file that does not have proper formatting.


# Contributions
Christina Pfeiffer and Gerlinde Grabmann contributed equally. 
ChatGPT and Copilot are acknowledged for code snipplets and generation of text used in Docstrings and the README file.

# Functionality implemented in the script
All functions of the assignment were implemented into the python script and extensively tested with the files provided. 

# Answer to the question: 

Using the files extension to check the file format can cause problems, when a file is given a wrong file extension by the user (e.g. FASTQ file with a .fasta extension). 
In that example the FASTQ file would be mistakenly treated as a fasta file leading to wrong results.

Our current approach is independent of the file extension and checks the file content itself. 
This could also lead to problems if there is a any other line before the sequences. This needs to be considered when using the script. 