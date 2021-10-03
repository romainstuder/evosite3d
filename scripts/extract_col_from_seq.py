#!/usr/bin/env python3

""" Extract transition between branches """

import sys
from Bio import SeqIO

# Load sequences
matrix = {}
# matrix_id = []

sequence_target = sys.argv[3]
seq_start = int(sys.argv[4])
seq_stop = int(sys.argv[5])
aln_start = 0
aln_stop = 0

max_length = 0
input_handle = open(sys.argv[1], 'r')
for record in SeqIO.parse(input_handle, 'fasta'):
    matrix[record.id] = str(record.seq)
    # matrix_id.append(record.id)
    if record.id == sequence_target:
        print(record.id)
        j = 0
        for i, aa in enumerate(str(record.seq)):
            if aa != "-":
                j = j+1
            if j == seq_start:
                aln_start = i-1
            elif j == seq_stop:
                aln_stop = i+1
                break

# Write file
input_handle = open(sys.argv[1], 'r')
file_out = open(sys.argv[2], "w")
for record in SeqIO.parse(input_handle, 'fasta'):
    matrix[record.id] = str(record.seq)
    file_out.write(">"+record.id+"\n")
    file_out.write(str(record.seq)[aln_start:aln_stop]+"\n")
file_out.close()
