#!/usr/bin/env python3

"""Extract probabilities"""

import argparse
from operator import itemgetter


def parse_probabilities(file_path, node):
    proba_dict = {}
    tag = False
    target_start = f"Prob distribution at node {node}, by site"
    target_end = f"Prob distribution at node {node + 1}, by site"

    with open(file_path, "r") as file_in:
        for line in file_in:
            line = line.rstrip()
            if target_end in line:
                tag = False
            if tag:
                tab = line.split()
                if len(tab) > 3:
                    site = tab[0]
                    probs = tab[3:24]
                    proba_dict[site] = probs
            if target_start in line:
                tag = True

    return proba_dict


def extract_site_probs(proba_dict):
    for site_str, prob_list in proba_dict.items():
        site = int(site_str)
        if 265 < site < 408:
            prob_aa_dict = {aa_prob[0]: float(aa_prob[2:7]) for aa_prob in prob_list}

            # Sort by probability descending
            sorted_probs = sorted(prob_aa_dict.items(), key=itemgetter(1), reverse=True)

            site_domain = site - 218
            aa, prob = sorted_probs[0]
            print(f"{site_domain} {aa}({prob:.2f})")


def main():
    parser = argparse.ArgumentParser(description="Extract probabilities from file at given node.")
    parser.add_argument("file_path", help="Path to the input file")
    parser.add_argument("node", type=int, help="Node number to extract probabilities for")
    args = parser.parse_args()

    proba_dict = parse_probabilities(args.file_path, args.node)
    extract_site_probs(proba_dict)


if __name__ == "__main__":
    main()
