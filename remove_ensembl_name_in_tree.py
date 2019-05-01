#!/usr/bin/env python3

"""Remove ensembl name from tree"""

import sys

def parsed_tree(tree):
    """

    Args:
        tree:

    Returns:

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


# Load tree
tree_line = ""
tree_file = open(sys.argv[1], "r")
while 1:
    line = tree_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tree_line = tree_line+ line
tree_file.close()

# Disassemble tree into list
tree_list = parsed_tree(tree_line)

# Write new tree
new_tree = ""
for item in tree_list:
    if "_" in item:
        item = item.split("_")[0]  # Keep only gene id
    new_tree = new_tree + item
new_tree = new_tree+";"  # Add final ";"
print(new_tree)
