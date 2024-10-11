"""Prokaryotic Gene Prediction Using Shine-Dalgarno Motif Detection.

Usage:
======
    python gpred/gpred.py -i [genome_file] -g [min_gene_len] -s [max_shine_dalgarno_distance] -d [min_gap] -p [predicted_genes_file] -o [fasta_output_file]

Arguments:
==========
    -i, --genome_file: (str) Complete genome file in fasta format.
    -g, --min_gene_len: (int) Minimum gene length to consider (default 50).
    -s, --max_shine_dalgarno_distance: (int) Maximum distance from start codon where to look for a Shine-Dalgarno motif (default 16).
    -d, --min_gap: (int) Minimum gap between two genes - shine box not included (default 40).
    -p, --predicted_genes_file: (str) Output file with the predicted genes.
    -o, --fasta_output_file: (str) Output file with the predicted genes in fasta format.

Example:
========
    python gpred/gpred.py -i data/listeria.fna -p results/predicted_genes_positions.csv -o results/predicted_genes.fasta

This command will predict genes in the Listeria genome file 'data/listeria.fna' with a minimum gene length of 50, a maximum distance of 16 between the start codon and the Shine-Dalgarno motif, a minimum gap of 40 between two genes, and will output the predicted genes in the files 'results/predicted_genes_positions.csv' and 'results/predicted_genes.fasta'.
"""

# LIBRARY IMPORTS
import argparse
import sys
import os
import csv
from typing import List, Union
from pathlib import Path
from textwrap import fill
import re
from re import Pattern
from loguru import logger


# METADATAS
__author__ = "Essmay Touami"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Essmay Touami"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Essmay Touami"
__email__ = "essmay.touami@etu.u-paris.fr"
__status__ = "Developpement"


# FUNCTIONS
def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments():  # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(
        description=__doc__, usage="{0} -h".format(sys.argv[0])
    )
    parser.add_argument(
        "-i",
        dest="genome_file",
        type=isfile,
        required=True,
        help="Complete genome file in fasta format",
    )
    parser.add_argument(
        "-g",
        dest="min_gene_len",
        type=int,
        default=50,
        help="Minimum gene length to consider (default 50).",
    )
    parser.add_argument(
        "-s",
        dest="max_shine_dalgarno_distance",
        type=int,
        default=16,
        help="Maximum distance from start codon "
        "where to look for a Shine-Dalgarno motif (default 16).",
    )
    parser.add_argument(
        "-d",
        dest="min_gap",
        type=int,
        default=40,
        help="Minimum gap between two genes - shine box not included (default 40).",
    )
    parser.add_argument(
        "-p",
        dest="predicted_genes_file",
        type=Path,
        default=Path("predict_genes.csv"),
        help="Tabular file giving position of predicted genes",
    )
    parser.add_argument(
        "-o",
        dest="fasta_file",
        type=Path,
        default=Path("genes.fna"),
        help="Fasta file giving sequence of predicted genes",
    )
    return parser.parse_args()


def read_fasta(fasta_file: Path) -> str:
    """Extract genome sequence from fasta files.

    :param fasta_file: (Path) Path to the fasta file.
    :return: (str) Sequence from the genome.
    """
    logger.info(f"Reading genome sequence from {fasta_file}...")
    with open(fasta_file, "r", encoding="utf-8") as file_in:
        sequence = []
        for line in file_in:
            if line.startswith(">"):
                continue
            sequence.append(
                line.strip().upper()
            )  # Ensure that the sequence is in uppercase

    seq = "".join(sequence)
    logger.debug(f"{len(seq)} nucleotides found.")
    logger.success("Genome sequence read successfully! \n")

    return seq


def find_start(
    start_regex: Pattern, sequence: str, start: int, stop: int
) -> Union[int, None]:
    """Find next start codon before a end position.

    :param start_regexp: A regex object that identifies a start codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :param stop: (int) Stop position of the research
    :return: (int) If exist, position of the start codon. Otherwise None.
    """
    # Search for a start codon in the sequence
    is_new_start_codon = start_regex.search(sequence, start, stop)
    if is_new_start_codon:
        # Return the position of the start codon
        return is_new_start_codon.start()
    return None


def find_stop(stop_regex: Pattern, sequence: str, start: int) -> Union[int, None]:
    """Find next stop codon that should be in the same reading phase as the start.

    :param stop_regexp: A regex object that identifies a stop codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :return: (int) If exist, position of the stop codon. Otherwise None.
    """
    # Search for a stop codon in the sequence after the start codon
    for i in range(start + 3, len(sequence) - 2, 3):
        codon = sequence[i : i + 3]
        if stop_regex.fullmatch(codon):
            # Return the position of the stop codon
            return i
    return None


def has_shine_dalgarno(
    shine_regex: Pattern, sequence: str, start: int, max_shine_dalgarno_distance: int
) -> bool:
    """Find a shine dalgarno motif before the start codon

    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Position of the start in the genome
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine dalgarno to the start position
    :return: (boolean) true -> has a shine dalgarno upstream to the gene, false -> no
    """
    # Define the search window
    pos_search = start - max_shine_dalgarno_distance
    # Check if the search window is valid
    if pos_search < 0:
        return False
    # Search for the motif in the search window
    for i in range(pos_search, start - 6):
        if shine_regex.search(sequence[i : i + len(shine_regex.pattern)]):
            return True

    return False


def predict_genes(
    sequence: str,
    start_regex: Pattern,
    stop_regex: Pattern,
    shine_regex: Pattern,
    min_gene_len: int,
    max_shine_dalgarno_distance: int,
    min_gap: int,
) -> List:
    """Predict most probable genes

    :param sequence: (str) Sequence from the genome.
    :param start_regexp: A regex object that identifies a start codon.
    :param stop_regexp: A regex object that identifies a stop codon.
    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param min_gene_len: (int) Minimum gene length.
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine dalgarno to the start position.
    :param min_gap: (int) Minimum distance between two genes.
    :return: (list) List of [start, stop] position of each predicted genes.
    """
    probable_genes = []
    current_pos = 0

    # Iterate over the sequence considering the minimum gap
    while len(sequence) - current_pos >= min_gap:
        # Find the next start and stop codon
        start_pos = find_start(start_regex, sequence, current_pos, len(sequence))
        stop_pos = find_stop(stop_regex, sequence, start_pos)
        # Checks if the start and stop codon are found
        if start_pos is None or stop_pos is None:
            current_pos += 1
            continue

        # Check if the gene is long enough
        if stop_pos - start_pos + 3 >= min_gene_len:
            # Check if there is a shine dalgarno motif before the start codon
            if has_shine_dalgarno(
                shine_regex, sequence, start_pos, max_shine_dalgarno_distance
            ):
                # Add the gene to the list
                probable_genes.append([start_pos + 1, stop_pos + 3])
                # Mve the current position after the gene added with the minimum gap
                current_pos = stop_pos + 3 + min_gap
            else:
                # Move to next codon if no shine dalgarno
                current_pos += 1
        else:
            # Move to next codon if gene is too short
            current_pos += 1

    logger.success(f"{len(probable_genes)} genes predicted successfully! \n")
    return probable_genes


def write_genes_pos(
    predicted_genes_file: Path, probable_genes: List[List[int]]
) -> None:
    """Write list of gene positions.

    :param predicted_genes_file: (Path) Output file of gene positions.
    :param probable_genes: List of [start, stop] position of each predicted genes.
    """
    logger.info(f"Writing predicted genes positions to {predicted_genes_file}...")

    # create results folder if it does not exist
    predicted_genes_file.parent.mkdir(parents=True, exist_ok=True)
    try:
        with predicted_genes_file.open("wt", encoding= "utf-8") as predict_genes:
            predict_genes_writer = csv.writer(predict_genes, delimiter=",")
            predict_genes_writer.writerow(["Start", "Stop"])
            predict_genes_writer.writerows(probable_genes)
        logger.success("Predicted genes positions written successfully! \n")
    except IOError:
        logger.error(f"Error cannot open {predicted_genes_file}")
        sys.exit(1)


def write_genes(
    fasta_file: Path,
    sequence: str,
    probable_genes: List[List[int]],
    sequence_rc: str,
    probable_genes_comp: List[List[int]],
):
    """Write gene sequence in fasta format

    :param fasta_file: (Path) Output fasta file.
    :param sequence: (str) Sequence of genome file in 5'->3'.
    :param probable_genes: (list) List of [start, stop] position of each predicted genes in 5'->3'.
    :param sequence_rc: (str) Sequence of genome file in 3' -> 5'.
    :param probable_genes_comp: (list)List of [start, stop] position of each predicted genes in 3' -> 5'.
    """
    logger.info(f"Writing predicted genes to {fasta_file}...")

    # create results folder if it does not exist
    fasta_file.parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(fasta_file, "wt") as fasta:
            for i, gene_pos in enumerate(probable_genes):
                fasta.write(
                    ">gene_{0}{1}{2}{1}".format(
                        i + 1, os.linesep, fill(sequence[gene_pos[0] - 1 : gene_pos[1]])
                    )
                )
            i = i + 1
            for j, gene_pos in enumerate(probable_genes_comp):
                fasta.write(
                    ">gene_{0}{1}{2}{1}".format(
                        i + 1 + j,
                        os.linesep,
                        fill(sequence_rc[gene_pos[0] - 1 : gene_pos[1]]),
                    )
                )
        logger.success("Predicted genes written successfully! \n")
    except IOError:
        logger.error(f"Error cannot open {fasta_file}")
        sys.exit(1)


def reverse_complement(sequence: str) -> str:
    """Get the reverse complement

    :param sequence: (str) DNA Sequence.
    :return: (str) Reverse complemented sequence.
    """
    logger.info(f"Reverse complementing the sequence ...{sequence[-10:]}")
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    reversed_seq = "".join([complement[base] for base in sequence[::-1]])

    logger.debug(f"Reversed sequence: {reversed_seq[:10]}...")
    logger.success("Reverse complement computed successfully! \n")
    return reversed_seq


# ==============================================================
# Main program
# ==============================================================
def main() -> None:  # pragma: no cover
    """
    Main program function
    """
    # Gene detection over genome involves to consider a thymine instead of
    # an uracile that we would find on the expressed RNA
    # start_codons = ['TTG', 'CTG', 'ATT', 'ATG', 'GTG']
    # stop_codons = ['TAA', 'TAG', 'TGA']
    start_regex = re.compile("AT[TG]|[ATCG]TG")
    stop_regex = re.compile("TA[GA]|TGA")
    # Shine AGGAGGUAA
    # AGGA ou GGAGG
    shine_regex = re.compile("A?G?GAGG|GGAG|GG.{1}GG")
    # Arguments
    args = get_arguments()

    # Read the genome sequence
    sequence = read_fasta(args.genome_file)

    # Predict genes
    logger.info(f"Predicting genes from {args.genome_file}...")
    probable_genes = predict_genes(
        sequence,
        start_regex,
        stop_regex,
        shine_regex,
        args.min_gene_len,
        args.max_shine_dalgarno_distance,
        args.min_gap,
    )

    # We reverse and complement the sequence to predict genes on the reverse strand
    sequence_rc = reverse_complement(sequence)
    logger.info("Predicting genes from reverse complemented sequence...")
    probable_genes_cp = predict_genes(
        sequence_rc,
        start_regex,
        stop_regex,
        shine_regex,
        args.min_gene_len,
        args.max_shine_dalgarno_distance,
        args.min_gap,
    )
    # Reverse the coordinates of the genes on the reverse strand
    probable_genes_cp = [
        [len(sequence_rc) - stop + 1, len(sequence_rc) - start + 1]
        for start, stop in probable_genes_cp
    ]

    # Write the predicted genes
    write_genes_pos(args.predicted_genes_file, probable_genes + probable_genes_cp)
    write_genes(
        args.fasta_file, sequence, probable_genes, sequence_rc, probable_genes_cp
    )


if __name__ == "__main__":
    main()
