#!/usr/bin/env python3

"""Extract full sequences from a formatted input file."""

import argparse


def extract_sequences(input_file):
    in_sequence_section = False

    with open(input_file, "r") as f:
        for line in f:
            line = line.rstrip()

            if "Overall accuracy of the 195 ancestral sequences" in line:
                in_sequence_section = False

            elif "List of extant and reconstructed sequences" in line:
                in_sequence_section = True
                continue

            if in_sequence_section and line and not line.startswith(" "):
                line = line.replace("e #", "e_")
                parts = line.split()
                if len(parts) > 3:
                    name = parts[0]
                    sequence = "".join(parts[1:])
                    print(f">{name}")
                    print(sequence)


def main():
    parser = argparse.ArgumentParser(
        description="Extract full sequences from a formatted input file."
    )
    parser.add_argument("input_file", help="Path to the input text file")
    args = parser.parse_args()

    extract_sequences(args.input_file)


if __name__ == "__main__":
    main()
