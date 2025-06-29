#!/usr/bin/env python3

"""Extract transition between branches in evolutionary sequence data."""

import argparse
import operator

from Bio import SeqIO


def parse_positions(seq, target_id, positions, shift):
    """Map alignment positions to sequence positions."""
    position_map = {}
    j = 0
    for i, aa in enumerate(seq):
        if j in positions:
            position_map[i] = j - shift
        if aa != "-":
            j += 1
    return position_map


def load_sequences(sequence_file, target_id, positions, shift):
    """Load sequences and find position mappings for the target sequence."""
    matrix = {}
    position_list = [int(pos) + shift for pos in positions.split("+")]
    position_map = {}

    with open(sequence_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            matrix[record.id] = str(record.seq)
            if record.id == target_id:
                position_map = parse_positions(str(record.seq), target_id, position_list, shift)

    return matrix, position_map


def extract_probabilities(rst_file, node, position_map):
    """Extract the top two most probable amino acids at specific positions."""
    site_list = []
    proba_list = []
    tag = False

    with open(rst_file, "r") as f:
        for line in f:
            line = line.strip()
            if f"Prob distribution at node {node + 1}, by site" in line:
                tag = False
            elif f"Prob distribution at node {node}, by site" in line:
                tag = True
            elif tag and line and len(line.split()) > 3 and not line.startswith("Prob"):
                tab = line.split()
                site = int(tab[0])
                if site in position_map:
                    probs = tab[3:24]
                    prob_dict = {entry[0]: float(entry[2:7]) for entry in probs}
                    top_probs = sorted(prob_dict.items(), key=operator.itemgetter(1), reverse=True)
                    site_list.append(position_map[site])
                    proba_list.append(top_probs[0])
                    print(position_map[site], top_probs[:2])

    return site_list, proba_list


def main():
    parser = argparse.ArgumentParser(
        description="Extract transition probabilities between branches in sequences."
    )
    parser.add_argument("sequence_file", help="FASTA file with sequences")
    parser.add_argument("rst_file", help="RST file with probability distributions")
    parser.add_argument("sequence_target", help="Target sequence ID to analyze")
    parser.add_argument("seq_of_interest", help="Positions of interest (format: pos1+pos2+pos3...)")
    parser.add_argument("shift", type=int, help="Shift to add to each position in seq_of_interest")
    parser.add_argument("node", type=int, help="Node number to extract probabilities for")

    args = parser.parse_args()

    _, position_map = load_sequences(
        args.sequence_file, args.sequence_target, args.seq_of_interest, args.shift
    )
    site_list, proba_list = extract_probabilities(args.rst_file, args.node, position_map)

    print("\t".join(map(str, site_list)))
    print("\t".join(f"{aa}({prob})" for aa, prob in proba_list))


if __name__ == "__main__":
    main()
