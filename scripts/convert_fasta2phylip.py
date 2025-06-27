#!/usr/bin/env python3

"""
Convert FASTA to PHYLIP format.
Usage: script.py input.fasta output.phy [name_length=50]
"""

import argparse
import sys
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description="Convert FASTA to PHYLIP format.")
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("output_phylip", help="Output PHYLIP file")
    parser.add_argument(
        "-n",
        "--name_length",
        type=int,
        default=50,
        help="Max length of sequence names (default: 50)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    print("Convert FASTA to PHYLIP")

    sequence_order = []
    sequences = {}

    with open(args.input_fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_id = record.id.split()[0]
            seq = str(record.seq).replace(" ", "")

            if "U" in seq:
                print(f"Sequence {seq_id} contains 'U'. Exiting.")
                sys.exit(1)

            sequence_order.append(seq_id)
            sequences[seq_id] = seq

    num_seqs = len(sequences)
    print(f"Number of sequences: {num_seqs}")

    lengths = {len(seq) for seq in sequences.values()}
    if len(lengths) != 1:
        print("Error: Sequences have differing lengths.")
        sys.exit(1)

    alignment_length = lengths.pop()
    print(f"Alignment length: {alignment_length}")
    print(f"Ratio (alignment length / 3): {alignment_length / 3:.2f}")

    if alignment_length % 3 != 0:
        print(
            "Warning: Alignment length is not divisible by 3; may not code for nucleotides."
        )

    with open(args.output_phylip, "w") as phyfile:
        phyfile.write(f"{num_seqs}\t{alignment_length}\n")
        for seq_id in sequence_order:
            name = args.name_length
            truncated_name = seq_id[:name].rstrip("_").replace(" ", "")
            phyfile.write(f"{truncated_name}  {sequences[seq_id]}\n")


if __name__ == "__main__":
    main()
