#!/usr/bin/env python3

"""Remove columns from FASTA alignment that are entirely gaps or 'X'."""

import argparse

from Bio import SeqIO


def remove_empty_columns(sequences):
    max_length = max(len(seq) for seq in sequences.values())

    filtered_seqs = {seq_id: [] for seq_id in sequences}

    for pos in range(max_length):
        column = [sequences[seq_id][pos] for seq_id in sequences]
        unique_chars = set(column)

        if len(unique_chars) == 1 and unique_chars.pop() in {"-", "X"}:
            continue

        for seq_id, char in zip(sequences, column):
            filtered_seqs[seq_id].append(char)

    return filtered_seqs


def parse_args():
    parser = argparse.ArgumentParser(
        description="Remove empty columns (all '-' or 'X') from a FASTA alignment."
    )
    parser.add_argument("input_fasta", help="Input FASTA alignment file")
    parser.add_argument("output_fasta", help="Output FASTA file without empty columns")
    return parser.parse_args()


def main():
    args = parse_args()
    sequences = {
        record.id: str(record.seq) for record in SeqIO.parse(args.input_fasta, "fasta")
    }
    filtered_seqs = remove_empty_columns(sequences)
    with open(args.output_fasta, "w") as out_handle:
        for seq_id, seq_chars in filtered_seqs.items():
            out_handle.write(f">{seq_id}\n{''.join(seq_chars)}\n")


if __name__ == "__main__":
    main()
