#!/usr/bin/env python3

"""
Convert FASTA to PHYLIP format.
Usage: script.py input.fasta output.phy [name_length=50]
"""

import argparse
import sys
from pathlib import Path

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


def validate_input_file(input_file):
    file_path = Path(input_file)
    if not file_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    if not file_path.is_file():
        raise ValueError(f"Path is not a file: {input_file}")


def load_sequences(input_fasta):
    sequence_order = []
    sequences = {}

    with open(input_fasta, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_id = record.id.split()[0]
            seq = str(record.seq).replace(" ", "")

            if "U" in seq:
                raise ValueError(f"Sequence {seq_id} contains 'U' (RNA nucleotide)")

            sequence_order.append(seq_id)
            sequences[seq_id] = seq

    if not sequences:
        raise ValueError("No sequences found in the input file")

    return sequence_order, sequences


def validate_alignment(sequences):
    lengths = {len(seq) for seq in sequences.values()}
    if len(lengths) != 1:
        raise ValueError("Sequences have differing lengths - not a proper alignment")

    return lengths.pop()


def write_phylip_file(output_phylip, sequence_order, sequences, alignment_length, name_length):
    num_seqs = len(sequences)

    with open(output_phylip, "w") as phyfile:
        phyfile.write(f"{num_seqs}\t{alignment_length}\n")
        for seq_id in sequence_order:
            truncated_name = seq_id[:name_length].rstrip("_").replace(" ", "")
            phyfile.write(f"{truncated_name}  {sequences[seq_id]}\n")


def main():
    try:
        args = parse_args()

        print("Convert FASTA to PHYLIP")

        validate_input_file(args.input_fasta)
        sequence_order, sequences = load_sequences(args.input_fasta)

        num_seqs = len(sequences)
        print(f"Number of sequences: {num_seqs}")

        alignment_length = validate_alignment(sequences)
        print(f"Alignment length: {alignment_length}")
        print(f"Ratio (alignment length / 3): {alignment_length / 3:.2f}")

        if alignment_length % 3 != 0:
            print("Warning: Alignment length is not divisible by 3; may not code for nucleotides.")

        write_phylip_file(
            args.output_phylip, sequence_order, sequences, alignment_length, args.name_length
        )
        print(f"Successfully converted to PHYLIP format: {args.output_phylip}")

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
