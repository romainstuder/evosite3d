#!/usr/bin/env python3
""" Extract transition between branches """

import sys

from Bio import SeqIO

# Load sequences
matrix = {}
# matrix_id = []

sequence_target = sys.argv[2]
seq_of_interest = sys.argv[3]
shift = int(sys.argv[4])
aln_start = 0
aln_stop = 0

position_str_list = seq_of_interest.split("+")
position_list = [int(x)+shift for x in position_str_list]


max_length = 0
input_handle = open(sys.argv[1], 'r')
for record in SeqIO.parse(input_handle, 'fasta'):
    matrix[record.id] = str(record.seq)
    # matrix_id.append(record.id)
    if record.id == sequence_target:
        j = 0
        for i, aa in enumerate(str(record.seq)):
            if j in position_list:
                print(record.id, j-shift, j, i)
            if aa != "-":
                j = j+1
