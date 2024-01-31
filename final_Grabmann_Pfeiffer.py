#! /usr/bin/env python3

import random
import gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO
import sys
import argparse
import logging
import shutil
import os

logging.basicConfig(level=logging.NOTSET)  # configure root logger
logger = logging.getLogger(__name__)  # create custom logger
# Logging levels: DEBUG/INFO/WARNING/ERROR/CRITICAL
logger.setLevel(logging.INFO)


def open_file(input_file):
    """Opens a file, which can be a gzipped file or a normal file.

    Args:
        input_file (str): Path to the file, which can be a gzipped file or a normal file.

    Returns:
        file handle: A handle to the opened file.

    Raises:
        ValueError: If the file does not exist or is not readable.
    """
    try:
        if input_file.endswith(".gz"):
            file = gzip.open(
                input_file, "rt"
            )  # rt = read text when handling compressed files
        else:
            file = open(input_file, "r")  # r = read
        # logger.info(f"Opened file: {input_file}")
        return file
    except Exception as e:
        logger.error(f"Error opening file {input_file}: {str(e)}")
        raise ValueError("The file does not exist or is not readable")


def detect_file_format(input_file):
    """ "Determine if the file is in FASTA, FASTQ format, or if the format is unknown.

    Args:
        input_file (str): Path to the file, which can be a gzipped file or a normal file.

    Returns:
        str: One of the following strings - 'fasta', 'fastq'.

    Raises:
        ValueError: If the file format is unknown
    """
    with open_file(input_file) as file:
        first_line = file.readline()
        if first_line.startswith(">"):
            file_format = "fasta"
        elif first_line.startswith("@"):
            file_format = "fastq"
        else:
            raise ValueError(f"Unknown file format {input_file}")
    # logger.info(f"Detected file format for {input_file}: {file_format}")
    return file_format


def read_records(input_file, use_fast_parsing=False):
    """Read sequences from file using SeqIO.index or fast parsing.

    Args:
        input_file (str): File to be handled.
        use_fast_parsing (bool): Decide if fast parsing should be used. Defaults to False.

    Returns:
        list or dict or SeqIO.index: If `use_fast_parsing` is True, returns a list of tuples
                                      (for FASTA) or a dictionary (for FASTQ) with the sequences.
                                      If False, returns a SeqIO.index object representing the
                                      sequences with identifiers as keys and corresponding
                                      sequences as values.
    """
    read_file = open_file(input_file)
    file_format = detect_file_format(input_file)

    records = None
    if use_fast_parsing:
        if file_format == "fasta":
            records = list(SimpleFastaParser(read_file))
        elif file_format == "fastq":
            records = list(FastqGeneralIterator(read_file))
    else:
        # Use SeqIO.index to read sequences from fasta file in case of large files
        # preventing to read all sequences into memory

        if input_file.endswith(
            ".gz"
        ):  # SeqIO.index cannot handle gzipped files, so decompress first
            file_name = decrompress_file(input_file)
            records = SeqIO.index(
                file_name, file_format
            )  # SeqIO.index returns a dictionary
        else:
            records = SeqIO.index(input_file, file_format)
            # SeqIO.index cannot deal with file handles, but needs a file name

    read_file.close()
    return records


def decrompress_file(input_file):
    """Decompresses a gzipped file and returns the file name. The original file is copied to a new file.
       This file and file name is needed for the SeqIO.index function.

    Args:
        input_file (str): Path to the file, which is a gzipped file.

    Returns:
        str: File name of the decompressed file as a string.

    Raises:
        ValueError: If an error occurs during decompression or copying, or if the file is not readable.
    """
    logger.info(f"Decompressing file: {input_file}")

    try:
        # Check if the destination file already exists
        output_file = input_file[:-3]
        if os.path.exists(output_file):
            user_input = input(
                f"Destination file '{output_file}' already exists. Do you want to overwrite it? (y/n): "
            ).lower()

            if user_input != "y":
                logger.warning(
                    f"Operation aborted. File '{output_file}' will not be overwritten."
                )
                print("Exiting script.")
                sys.exit()
        # Copy the original file to a new file
        with gzip.open(input_file, "rt") as input_f:
            with open(output_file, "w") as output_f:
                shutil.copyfileobj(input_f, output_f)

        file_name = output_file

    except Exception as e:
        logger.error(f"Error decompressing file {input_file}: {str(e)}")
        raise ValueError("The file does not exist or is not readable")

    return file_name


def subsample_sequences(
    input_file, num_reads, use_fast_parsing=False, set_seed=False, preserve_format=False
):
    """Subsample sequences from a file.

    Args:
        input_file (str): Path to the file containing sequences.
        num_reads (int): Number of reads to subsample.
        use_fast_parsing (bool, optional): Decide whether to use fast parsing methods. Defaults to False.
        set_seed (bool, optional): If True, set a seed for reproducibility. Defaults to False.
        preserve_format (bool, optional): If True, preserve the original formatting of the sequences. Defaults to False.

    Raises:
        ValueError: If the number of reads is larger than the number of sequences in the file.

    Returns:
        None. Writes to stdout or calls fasta_output function.
    """

    if set_seed:
        random.seed(set_seed)

    records = read_records(input_file, use_fast_parsing)

    num_sequences = len(records)
    assert (
        num_sequences >= num_reads
    ), f"Cannot subsample {num_reads} sequences from file {input_file} with {num_sequences} sequences"

    if use_fast_parsing:
        if not preserve_format:
            selected_records = random.sample(records, num_reads)
            fasta_output(input_file, selected_records, use_fast_parsing)

    else:
        selected_keys = random.sample(list(records.keys()), num_reads)
        selected_records = {key: records[key] for key in selected_keys}
        if preserve_format:
            logger.info("Preserve format is needed; subsampling is printed to stdout")
            for seq_id in selected_records.keys():
                formatted_sequence = records.get_raw(seq_id).decode()
                sys.stdout.write(formatted_sequence)
        else:
            fasta_output(input_file, selected_records, use_fast_parsing)


def fasta_output(input_file, selected_records, use_fast_parsing):
    """Print subsampled sequences to stdout in FASTA or FASTQ format not preserving the number of lines.

    Args:
        input_file (str): Path to the original file containing sequences.
        selected_records (list or dict): Subsampled records to be printed. If fast parsing is used, a list of tuples is expected. If not, a dictionary is expected.
        use_fast_parsing (bool): Indicates whether fast parsing methods were used.




    Returns:
        None. Writes to stdout.
    """
    file_format = detect_file_format(input_file)
    if use_fast_parsing:
        logger.info("Preserve format is not needed; subsampling is printed to stdout")
        if file_format == "fasta":
            for title, seq in selected_records:
                sys.stdout.write(f">{title}\n{seq}\n")
        elif file_format == "fastq":
            for title, seq, qual in selected_records:
                sys.stdout.write(f"@{title}\n{seq}\n+\n{qual}\n")

    else:
        logger.info("Preserve format is not needed; subsampling is printed to stdout")
        for seq_id in selected_records.keys():
            record_seq_id = selected_records[seq_id]
            formatted_sequence = record_seq_id.format(file_format)
            sys.stdout.write(f"{formatted_sequence}")


def main(args):
    """Execute the main functionality of the subsampling program.

    Args:
        args (argparse.Namespace): Command-line arguments parsed by argparse.

    Returns:
        None: The function prints subsampled sequences to stdout.

    Usage:
        This function is intended to be called from the command line with arguments
        specifying the input file, the number of reads to subsample, and optional flags
        for fast parsing, preserving format, and setting a seed for random number generation.

    Example:
        $ python final_Grabmann_Pfeiffer.py fasta_testfile.fasta 3 -f -p -s 42

    Notes:
        - If both -f (fast) and -p (preserve format) flags are used together, the combination is not allowed.
          Please use the -p flag without the -f flag.

    """

    logger.debug("Executing main function.")

    if args.preserve_format and args.fast:
        logger.info(
            "Preserve format is true and fast is true: \n This combination is not possible and cannot be executed. \n Please use the -p flag without the -f flag."
        )
        sys.exit(1)

    subsample_sequences(
        args.input_file, args.num_reads, args.fast, args.seed, args.preserve_format
    )


if __name__ == "__main__":
    # Parse and check arguments
    parser = argparse.ArgumentParser(
        description="""Parse (compressed) fasta/q files and and randomly select sequences to generate a
        subsample ."""
    )

    parser.add_argument(
        "input_file",
        metavar="FILE",
        action="store",
        type=str,
        help="""Input file (fasta or fastq, can be gzipped)""",
    )

    parser.add_argument(
        "num_reads",
        metavar="NUM",
        action="store",
        type=int,
        help="""Number of reads to subsample""",
    )

    parser.add_argument(
        "-f",
        "--fast",
        action="store_true",
        help="""Use fast parsing (only for fasta/fastq, not for gzipped files)""",
    )

    parser.add_argument(
        "-p",
        "--preserve-format",
        action="store_true",
        help="""Preserve original formatting of sequences""",
    )

    parser.add_argument(
        "-s",
        "--seed",
        metavar="NUM",
        action="store",
        type=int,
        help="""Seed for random number generator""",
    )

    args = parser.parse_args()

    # Check if number of reads is positive
    if args.num_reads <= 0:
        logger.error(f"Number of reads must be positive.")
        sys.exit(1)

    # logger.warning(f"Program arguments: {args}")
    main(args)
