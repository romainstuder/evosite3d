#!/usr/bin/env python3

"""Get position(s) in CDS from Trimal output for multiple site indices."""

import argparse

from Bio import SeqIO


def load_sequences(fasta_path):
    sequences = {}
    with open(fasta_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            protein_seq = record.seq.translate(to_stop=True, gap="-")
            sequences[record.id] = protein_seq
    return sequences


def load_trimal_columns(trimal_file_path):
    column_str_list = []
    with open(trimal_file_path, "r") as file_in:
        for line in file_in:
            line = line.rstrip()
            if line.startswith("#"):
                cleaned_line = line.replace(",", "")
                column_str_list = cleaned_line.split()[1:]
    return [int(x) for x in column_str_list]


def build_position_mapping(sequences):
    position_mappings = {}
    for gene, seq in sequences.items():
        mapping = {}
        cds_pos = 0
        for align_pos, aa in enumerate(seq, start=1):
            if aa != "-":
                cds_pos += 1
                mapping[align_pos] = cds_pos
        position_mappings[gene] = mapping
    return position_mappings


def parse_args():
    parser = argparse.ArgumentParser(description="Map Trimal output columns to CDS positions.")
    parser.add_argument("full_alignment", help="Path to the full (untrimmed) alignment FASTA file")
    parser.add_argument("trimal_output", help="Path to the Trimal output file with kept columns")
    parser.add_argument("reference_gene", help="Reference gene ID present in the alignment")
    parser.add_argument(
        "site_indices",
        type=str,
        help="Site indices (0-based) separated by comma or space",
    )
    return parser.parse_args()


def parse_site_indices(site_indices_raw):
    # Split by comma or whitespace
    if "," in site_indices_raw:
        tokens = site_indices_raw.split(",")
    else:
        tokens = site_indices_raw.split()
    return [int(x) for x in tokens]


def main():
    args = parse_args()

    sequences = load_sequences(args.full_alignment)
    column_list = load_trimal_columns(args.trimal_output)
    position_mappings = build_position_mapping(sequences)

    if args.reference_gene not in sequences:
        raise KeyError(f"Reference gene '{args.reference_gene}' not found in sequences.")

    site_indices = parse_site_indices(args.site_indices)
    max_site_index = len(column_list) // 3 - 1

    for site_index in site_indices:
        if not (0 <= site_index <= max_site_index):
            print(f"Warning: site_index {site_index} out of range (0 to {max_site_index})")
            continue

        raw_col = column_list[site_index * 3]
        column_in_full = raw_col // 3

        cds_pos = position_mappings[args.reference_gene].get(column_in_full)
        if cds_pos is None:
            print(
                f"Warning: No CDS position found for column {column_in_full} "
                f"in gene {args.reference_gene}"
            )
            cds_pos = "N/A"

        aa = sequences[args.reference_gene][column_in_full - 1] if column_in_full > 0 else "N/A"
        print(f"site_index: {site_index}\tCDS pos: {cds_pos}\tAA: {aa}")


if __name__ == "__main__":
    main()
