# Subsample_biological_sequences

Subsampling Tool

Overview

This Python script is designed to parse (compressed) FASTA or FASTQ files and randomly select sequences to generate a subsample.
It provides options for fast parsing, preserving the original formatting of sequences, and setting a seed for reproducibility.

Usage

Requirements
    Python 3.9
    Biopython library

Installation
    Clone the repository:
    git clone https://github.com/your_username/your_repo.git

Install the required dependencies:
    pip install biopython ##macht man das so??

Running the Script
Use the script from the command line with the following format:

python Subsample_CP_LG.py INPUT_FILE NUM_READS [-f] [-p] [-s SEED] ####name

INPUT_FILE: Path to the input file (FASTA or FASTQ, can be gzipped).
NUM_READS: Number of reads to subsample.
-f or --fast: Use fast parsing.
-p or --preserve-format: Preserve the original formatting of sequences.
-s SEED or --seed SEED: Seed for the random number generator for reproducibility or paired-end reads.

Examples
# Subsample 3 reads from a FASTA file with fast parsing, preserving format, and setting seed
python Subsample_CP_LG_edit.py fasta_testfile.fasta 3 -f -p -s 42

# Subsample 3 reads from a FASTQ file with fast parsing, preserving format, and setting seed
python Subsample_CP_LG_edit.py Test_L001_R1_001.fastq 3 -f -p -s 42

# Subsample 3 reads from a gzipped FASTQ file with fast parsing, preserving format, and setting seed
python Subsample_CP_LG_edit.py Test_L001_R1_001.fastq.gz 3 -f -p -s 42

Notes

The script supports both FASTA and FASTQ file formats, and it can handle gzipped files.
All functions can be called as stand alone modules.

Large files that do not fit into memory should not use the fast option. 
Original format cannot be preserved for fast parsing.

The -f flag enables fast parsing for large files that fit into memory.
The -p flag preserves the original formatting of sequences.
The -s flag allows setting a seed for reproducibility.

Contributions
Christina Pfeiffer and Gerlinde Grabmann contributed equally. 
ChatGPT and Copilot are acknowledged for code snippets and phrasing. 