#!/usr/bin/env python3

""" Remove empty columns"""

import sys
from Bio import SeqIO

# Load sequences
matrix = {}
# matrix_id = []

max_length = 0
input_handle = open(sys.argv[1], 'r')
for record in SeqIO.parse(input_handle, 'fasta'):
    matrix[record.id] = str(record.seq)
    # matrix_id.append(record.id)
    if len(str(record.seq)) > max_length:
        max_length = len(str(record.seq))

matrix_new = {}
for gene_id in matrix:
    matrix_new[gene_id] = ""

for pos in range(max_length):
    aa_col = []
    for gene_id in matrix:
        aa = matrix[gene_id][pos]
        aa_col.append(aa)

    if not ((len((set(aa_col))) == 1) and (list(set(aa_col))[0] in ["-", "X"])):
        for gene_id in matrix:
            seq = matrix_new[gene_id]
            seq = seq+matrix[gene_id][pos]
            matrix_new[gene_id] = seq

# Write file
file_out = open(sys.argv[2], "w")
for gene_id in matrix:
    file_out.write(">"+gene_id+"\n")
    file_out.write("".join(matrix_new[gene_id])+"\n")
file_out.close()
