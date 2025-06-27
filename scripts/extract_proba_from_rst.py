#!/usr/bin/env python3

"""Extract transition between branches"""

import sys
import operator

from Bio import SeqIO

# Load sequences
matrix = {}
# matrix_id = []

sequence_file = sys.argv[1]
rst_file = sys.argv[2]
sequence_target = sys.argv[3]
seq_of_interest = sys.argv[4]
shift = int(sys.argv[5])
node = int(sys.argv[6])
aln_start = 0
aln_stop = 0

position_str_list = seq_of_interest.split("+")
position_list = [int(x) + shift for x in position_str_list]


list_of_position = {}

max_length = 0
input_handle = open(sequence_file, "r")
for record in SeqIO.parse(input_handle, "fasta"):
    matrix[record.id] = str(record.seq)
    # matrix_id.append(record.id)
    if record.id == sequence_target:
        j = 0
        for i, aa in enumerate(str(record.seq)):
            if j in position_list:
                # print(record.id, j-shift, j, i)
                list_of_position[i] = j - shift
            if aa != "-":
                j = j + 1


# Load rst

tag = 0
branch = ""
file_in = open(rst_file, "r")

site_list = []
proba_list = []

while 1:
    line = file_in.readline()
    if line == "":
        break
    line = line.rstrip()
    if "Prob distribution at node " + str(node + 1) + ", by site" in line:
        tag = 0
    if tag == 1:
        tab = line.split()
        if len(tab) > 3:
            # print(tab)
            if tab[0] != "Prob":
                site = int(tab[0])
                if site in list_of_position:
                    prob_list = tab[3:24]
                    prob_dict = {}
                    for aa_prob in prob_list:
                        prob_dict[aa_prob[0]] = float(aa_prob[2:7])

                    sorted_x = sorted(
                        prob_dict.items(), key=operator.itemgetter(1), reverse=True
                    )

                    print(list_of_position[site], sorted_x[0:2])
                    site_list.append(list_of_position[site])
                    proba_list.append(sorted_x[0])
    if "Prob distribution at node " + str(node) + ", by site" in line:
        tag = 1

file_in.close()

site_str_list = [str(x) for x in site_list]
proba_str_list = [str(x[0]) + "(" + str(x[1]) + ")" for x in proba_list]

print("\t".join(site_str_list))
print("\t".join(proba_str_list))
