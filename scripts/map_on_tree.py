#!/usr/bin/env python3
"""Module to map info on tree"""

import argparse
from typing import List

from Bio import SeqIO
from Bio.SeqUtils import ProtParam


def load_ancestral_sequences(file_path: str) -> dict:
    """Parse a FASTA file and return a dictionary of node IDs to isoelectric points."""
    node_dict = {}
    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = str(record.seq).replace("-", "")
            analysed_protein = ProtParam.ProteinAnalysis(sequence)
            pI = analysed_protein.isoelectric_point()
            # MW = analysed_protein.molecular_weight()
            node_dict[record.id] = round(pI, 2)
    return node_dict


def load_tree(file_path: str) -> List[str]:
    """Read and tokenise the tree file."""
    with open(file_path, "r") as tree_file:
        return tree_file.read().split()


def map_iso_point_to_tree(tree_tab: List[str], node_dict: dict) -> str:
    new_tree_line = ""
    for item in tree_tab:
        if "node" + item in node_dict:
            # print "node"+item, node_dict["node"+item]
            iso_point = node_dict["node" + item]
            item = str(round(iso_point, 2))
        new_tree_line = new_tree_line + item
    return new_tree_line


# def map_pI_to_tree(tree_tab: List[str], node_dict: dict) -> str:
#     """Replace node tokens in the tree with their corresponding pI values."""
#     return "".join(
#         str(node_dict.get(f"node{token}", token)) for token in tree_tab
#     )


def main():
    parser = argparse.ArgumentParser(description="Map isoelectric points to tree nodes")
    parser.add_argument("ancestral_fasta", help="Path to ancestral FASTA file")
    parser.add_argument("tree_file", help="Path to tree file")

    args = parser.parse_args()

    node_dict = load_ancestral_sequences(args.ancestral_fasta)
    tree_tokens = load_tree(args.tree_file)
    new_tree = map_iso_point_to_tree(tree_tokens, node_dict)

    print(new_tree)


if __name__ == "__main__":
    main()
