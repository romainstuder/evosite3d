#!/usr/bin/env python3

"""Translate DNA sequences from a FASTA file to protein sequences."""

import argparse

from Bio import SeqIO


def translate_fasta(file_path):
    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            print(f">{record.id}")
            print(record.seq.translate(to_stop=True))


def main():
    parser = argparse.ArgumentParser(
        description="Translate DNA sequences from a FASTA file to protein sequences."
    )
    parser.add_argument("fasta_file", help="Path to the input FASTA file containing DNA sequences")
    args = parser.parse_args()

    translate_fasta(args.fasta_file)


if __name__ == "__main__":
    main()
