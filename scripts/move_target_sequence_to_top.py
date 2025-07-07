#!/usr/bin/env python3

import argparse

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


def move_target_sequence_to_top(input_file, output_file, target_id, file_format="fasta"):
    # Load the alignment
    alignment = AlignIO.read(input_file, file_format)

    # Find the target sequence
    target_record = None
    other_records = []

    for record in alignment:
        if target_id in record.id:
            target_record = record
        else:
            other_records.append(record)

    if target_record is None:
        raise ValueError(f"Target ID '{target_id}' not found in the alignment.")

    # Reorder the alignment
    new_alignment = MultipleSeqAlignment([target_record] + other_records)

    # Write the output
    AlignIO.write(new_alignment, output_file, file_format)
    print(f"Saved reordered alignment to '{output_file}' with '{target_id}' on top.")


def main():
    parser = argparse.ArgumentParser(description="Move a target sequence to the top of a MSA file.")
    parser.add_argument("input_file", help="Input MSA file (e.g., alignment.fasta)")
    parser.add_argument("output_file", help="Output MSA file")
    parser.add_argument("target_id", help="Target sequence ID or substring to match")
    parser.add_argument("-f", "--format", default="fasta", help="MSA file format (default: fasta)")

    args = parser.parse_args()
    move_target_sequence_to_top(args.input_file, args.output_file, args.target_id, args.format)


if __name__ == "__main__":
    main()
