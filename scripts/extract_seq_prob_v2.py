#!/usr/bin/env python3
"""Extract transition between branches"""

import operator
import sys

from Bio import SeqIO

# Load sequences
matrix = {}
# matrix_id = []

sequence_file = sys.argv[1]
rst_file = sys.argv[2]
sequence_target = sys.argv[3]
# seq_of_interest = sys.argv[4]
shift = int(sys.argv[5])
node = int(sys.argv[6])
aln_start = 0
aln_stop = 0

# position_list = seq_of_interest.split("+")
# position_list = [int(x)+shift for x in position_list]


list_of_position = {}

sites_found = []

max_length = 0
input_handle = open(sequence_file, "r")
for record in SeqIO.parse(input_handle, "fasta"):
    matrix[record.id] = str(record.seq)
    # matrix_id.append(record.id)
    if record.id == sequence_target:
        j = 0
        for i, aa in enumerate(str(record.seq)):
            if (
                (j - shift > 48)
                and (j - shift < 179)
                and (j - shift not in sites_found)
            ):
                sites_found.append(j - shift)
                # print(record.id, j-shift, j, i)
                list_of_position[i] = j - shift
                if j - shift == 178:
                    break
            if aa != "-":
                j = j + 1

# sys.exit()
# print(list_of_position)
# Load rst

tag = 0
branch = ""
file_in = open(rst_file, "r")


site_list_1 = []
prob_list_1 = []


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

                    # print(list_of_position[site], site,
                    #       str(sorted_x[0][0]) +
                    #       "(" + str(round(sorted_x[0][1], 2)).ljust(4, "0") + ")")
                    print(
                        list_of_position[site],
                        str(sorted_x[0][0])
                        + "("
                        + str(round(sorted_x[0][1], 2)).ljust(4, "0")
                        + ")",
                    )
                    site_list_1.append(str(list_of_position[site]))
                    prob_list_1.append(
                        str(sorted_x[0][0])
                        + "("
                        + str(round(sorted_x[0][1], 2)).ljust(4, "0")
                        + ")"
                    )
                    # sorted_x[0:3])
                    # print(str(site_domain),
                    #       str(sorted_x[0][0])+"("+str(round(sorted_x[0][1], 2)).ljust(4, "0")+")")

    if "Prob distribution at node " + str(node) + ", by site" in line:
        tag = 1

file_in.close()

print(len(site_list_1))
print(len(prob_list_1))

for i in range(len(site_list_1)):
    print(site_list_1[i])
    print(prob_list_1[i])

print("\t".join(site_list_1))
print("\t".join(prob_list_1))

file_out = open("tmp.txt", "w")
file_out.write("\t".join(site_list_1) + "\n")
file_out.write("\t".join(prob_list_1) + "\n")
file_out.close()
