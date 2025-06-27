#!/usr/bin/env python3

"""Extract sequences from a FASTA file based on a list of gene IDs."""

import sys
from Bio import SeqIO


def load_gene_list(filepath):
    """Read a list of gene IDs from a file."""
    with open(filepath, "r") as f:
        return {line.strip().split()[0] for line in f if line.strip()}


def extract_sequences(fasta_file, gene_set):
    """Yield sequences from a FASTA file whose IDs are in the gene_set."""
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in gene_set:
                yield record


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <gene_list.txt> <sequences.fasta>")
        sys.exit(1)

    gene_list_file = sys.argv[1]
    fasta_file = sys.argv[2]

    gene_set = load_gene_list(gene_list_file)
    for record in extract_sequences(fasta_file, gene_set):
        print(f">{record.id}")
        print(str(record.seq))


if __name__ == "__main__":
    main()
