#!/usr/bin/env python3

"""Compute the isoelectric point (pI) and molecular weight (MW) of proteins from a FASTA file."""

# Load libraries
import argparse
from Bio import SeqIO
from Bio.SeqUtils import ProtParam


def compute_properties(fasta_file):
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = str(record.seq).replace("-", "")
            analysed_protein = ProtParam.ProteinAnalysis(sequence)

            pI = analysed_protein.isoelectric_point()
            MW = analysed_protein.molecular_weight()

            print(f"{record.id}\t{pI:.2f}\t{MW:.2f}")


def main():
    parser = argparse.ArgumentParser(description="Compute pI and MW from a FASTA file.")
    parser.add_argument("fasta_file", help="Input FASTA file with protein sequences.")
    args = parser.parse_args()

    compute_properties(args.fasta_file)


if __name__ == "__main__":
    main()
