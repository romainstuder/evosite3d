#!/usr/bin/env python3

"""Extract transition between branches"""

import argparse

# target_branches = ["197..199", "199.200",
#                    "208..209",
#                    "209..329",
#                    "209..210",
#                    "210..265",
#                    "210..211", "211..212", "212..213"]

target_branches = ["208..209", "209..210", "209..265", "208..230"]


# target_position = ["239", "245", "272", "282",
#                    "299", "318", "320", "355",
#                    "366"]

target_position = ["198", "204", "222", "234", "251", "270", "272", "307", "318"]


def extract_transitions(file_path):
    tag = 0
    branch = ""
    with open(file_path, "r") as file_in:
        for line in file_in:
            line = line.rstrip()
            if "List of extant and reconstructed sequences" in line:
                tag = 0
            if tag == 1:
                tab = line.split()
                if tab:
                    if tab[0] == "Branch":
                        branch = tab[2]
                    position = tab[0]
                    if position in target_position and branch in target_branches:
                        print(branch, position, line)
            if "Summary of changes along branches" in line:
                tag = 1


def main():
    parser = argparse.ArgumentParser(description="Extract transitions between branches from a file")
    parser.add_argument("input_file", help="Input file to parse")
    args = parser.parse_args()

    extract_transitions(args.input_file)


if __name__ == "__main__":
    main()
