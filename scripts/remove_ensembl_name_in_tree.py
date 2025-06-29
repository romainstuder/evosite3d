#!/usr/bin/env python3

"""Remove ensembl name from tree"""

import argparse


def parsed_tree(tree):
    """
    Parse a Newick-formatted tree into components while extracting parts.

    Args:
        tree (str): The raw Newick tree string.

    Returns:
        list: Components of the tree split into a list.
    """
    tab_tree = []
    word = ""
    tag = 0
    tag_bracket = 0
    for char in tree:
        if tag == 0:
            if char == "(":
                tag = 1
        if char == "[":
            if word:
                tab_tree.append(word)
            word = char
            tag_bracket = 1
        elif char == "]":
            word = word + char
            tab_tree.append(word)
            word = ""
            tag_bracket = 0
        elif tag_bracket == 1:
            word = word + char
        elif tag == 1:
            if char in ("(", ")", ","):
                tag = 1
                if word != "":
                    tab_tree.append(word)
                word = char
                tab_tree.append(word)
                word = ""
            elif char == ":":
                if word != "":
                    tab_tree.append(word)
                word = char
            elif char == " ":
                tab_tree.append(word)
                word = ""
            elif (char == "." or "_") and (char != ";"):  # Any alphabetical character
                word = word + char
            elif char == ";":
                if word != "":
                    tab_tree.append(word)
                # print "End of tree reach"
                break
    return tab_tree


def clean_tree(input_path):
    """Read, parse, and clean the tree by removing Ensembl names.

    Args:
        input_path (str): Path to the file containing the tree.

    Returns:
        str: Cleaned Newick tree string.
    """
    with open(input_path, "r") as f:
        tree_line = "".join(line.strip() for line in f)

    # Disassemble tree into list
    tree_list = parsed_tree(tree_line)

    # Write new tree
    new_tree = ""
    for item in tree_list:
        if "_" in item:
            item = item.split("_")[0]  # Keep only gene id
        new_tree += item
    return new_tree + ";"  # Add final ";"


def main():
    parser = argparse.ArgumentParser(description="Remove Ensembl names from a Newick tree.")
    parser.add_argument("input_file", help="Path to the input Newick tree file")
    args = parser.parse_args()

    cleaned = clean_tree(args.input_file)
    print(cleaned)


if __name__ == "__main__":
    main()
