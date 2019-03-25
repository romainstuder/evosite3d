#!/usr/bin/env python3

"""Extract sequences"""

import sys
from Bio import SeqIO

# Load list of sequences to keep
sequence_list = []
file_in = open(sys.argv[1], "r")
while 1:
    line = file_in.readline()
    if line == "":
        break
    tab = line.split()
    gene = tab[0]
    sequence_list.append(gene)
file_in.close()

# Write only sequences that are in the list to keep
for record in SeqIO.parse(open(sys.argv[2], "rU"), "fasta"):
    gene = record.id
    if gene in sequence_list:
        sequence = str(record.seq)
        print(">"+gene)
        print(sequence)
