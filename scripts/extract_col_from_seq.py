#!/usr/bin/env python3

"""Extract a sequence region between specified positions in an alignment, accounting for gaps."""

import argparse
import sys

from Bio import SeqIO


def find_alignment_coordinates(sequence, seq_start, seq_stop):
    """Map ungapped sequence positions to alignment coordinates in the alignment."""
    j = 0
    aln_start = 0
    aln_stop = 0
    for i, aa in enumerate(sequence):
        if aa != "-":
            j += 1
        if j == seq_start:
            aln_start = i - 1
        elif j == seq_stop:
            aln_stop = i + 1
            break
    return aln_start, aln_stop


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Extract an alignment region based on ungapped sequence positions."
    )
    parser.add_argument("input_fasta", help="Input FASTA file containing alignment")
    parser.add_argument("output_fasta", help="Output FASTA file for extracted region")
    parser.add_argument("target_sequence_id", help="ID of the reference sequence in the alignment")
    parser.add_argument("start_pos", type=int, help="Start position (ungapped, 1-based)")
    parser.add_argument("stop_pos", type=int, help="Stop position (ungapped, 1-based)")
    return parser.parse_args()


def main():
    args = parse_arguments()

    records = list(SeqIO.parse(args.input_fasta, "fasta"))
    sequence_dict = {record.id: str(record.seq) for record in records}

    if args.target_sequence_id not in sequence_dict:
        sys.exit(f"Error: Sequence ID '{args.target_sequence_id}' not found in input.")

    target_sequence = sequence_dict[args.target_sequence_id]
    aln_start, aln_stop = find_alignment_coordinates(target_sequence, args.start_pos, args.stop_pos)

    if aln_start is None or aln_stop is None:
        sys.exit("Error: Could not map specified sequence positions to alignment.")

    with open(args.output_fasta, "w") as out_handle:
        for record in records:
            subseq = str(record.seq)[aln_start:aln_stop]
            out_handle.write(f">{record.id}\n{subseq}\n")


if __name__ == "__main__":
    main()
