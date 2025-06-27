#!/usr/bin/env python3

"""Translate DNA sequences from a FASTA file to protein sequences."""

import sys
from Bio import SeqIO


def translate_fasta(file_path):
    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            print(f">{record.id}")
            print(record.seq.translate(to_stop=True))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <fasta_file>", file=sys.stderr)
        sys.exit(1)
    translate_fasta(sys.argv[1])
