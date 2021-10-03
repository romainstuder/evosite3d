#!/usr/bin/env python3
""" Module to map info on tree """

import sys
from typing import List

from Bio import SeqIO
from Bio.SeqUtils import ProtParam

# Load ancestral sequences
node_dict = {}
handle = open(sys.argv[1], "rU")
for record in SeqIO.parse(handle, "fasta"):
    sequence = record.seq.tostring()
    sequence = sequence.replace("-", "")
    analysed_protein = ProtParam.ProteinAnalysis(sequence)

    # Compute some properties
    pI = analysed_protein.isoelectric_point()
    MW = analysed_protein.molecular_weight()
    # print record.id, pI, MW
    node_dict[record.id] = pI
handle.close()


# Load tree
tree_tab: List[str] = []
tree_file = open(sys.argv[2], "r")
while 1:
    line = tree_file.readline()
    if line == "":
        break
    tab = line.split()
    tree_tab = tree_tab+tab
tree_file.close()

# print tree_tab

new_tree_line = ""
for item in tree_tab:
    if "node"+item in node_dict:
        # print "node"+item, node_dict["node"+item]
        pI = node_dict["node"+item]
        pI = round(pI, 2)
        item = str(pI)
    new_tree_line = new_tree_line+item

print(new_tree_line)
