#!/usr/bin/env python3
"""Extract transition between branches"""

import argparse
from Bio import SeqIO
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract transition positions between branches from a FASTA alignment."
    )
    parser.add_argument(
        "fasta_file", type=str, help="Input FASTA file with aligned sequences."
    )
    parser.add_argument(
        "sequence_target", type=str, help="Target sequence ID to analyze."
    )
    parser.add_argument(
        "seq_of_interest",
        type=str,
        help="Positions of interest separated by '+', e.g., '10+20+30'.",
    )
    parser.add_argument("shift", type=int, help="Integer shift to apply to positions.")
    return parser.parse_args()


def main():
    args = parse_args()

    position_str_list = args.seq_of_interest.split("+")
    position_list = [int(pos) + args.shift for pos in position_str_list]

    matrix = {}

    with open(args.fasta_file, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            matrix[record.id] = str(record.seq)

    seq = matrix.get(args.sequence_target)
    if seq is None:
        print(
            f"Error: Sequence ID '{args.sequence_target}' not found in the input file.",
            file=sys.stderr,
        )
        sys.exit(1)

    j = 0
    for i, aa in enumerate(seq):
        if j in position_list:
            print(args.sequence_target, j - args.shift, j, i)
        if aa != "-":
            j += 1


if __name__ == "__main__":
    main()
