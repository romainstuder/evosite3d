#!/usr/bin/env python3

"""Extract and reformat node information from an RST file."""

import argparse


def parse_rst_file(input_file):
    with open(input_file, "r") as file_in:
        for line in file_in:
            line = line.rstrip()
            if line.startswith("node"):
                tab = line.split("          ")  # Deliberate spacing used as delimiter
                if len(tab) >= 2:
                    node_id = tab[0].replace(" ", "").replace("#", "")
                    sequence = tab[1].replace(" ", "")
                    print(f">{node_id}")
                    print(sequence)


def main():
    parser = argparse.ArgumentParser(
        description="Extract sequences from RST file nodes."
    )
    parser.add_argument("input_file", help="Path to the RST file.")
    args = parser.parse_args()
    parse_rst_file(args.input_file)


if __name__ == "__main__":
    main()
