#!/usr/bin/env python3

"""Extract probabilities"""

import operator
import sys

TAG = 0
BRANCH = ""

file_in = open(sys.argv[1], "r")
node = int(sys.argv[2])

proba_dict = {}

while 1:
    line = file_in.readline()
    if line == "":
        break
    line = line.rstrip()
    if "Prob distribution at node " + str(node + 1) + ", by site" in line:
        TAG = 0
    if TAG == 1:
        tab = line.split()
        if len(tab) > 3:
            # print(tab)
            site = tab[0]
            prob = tab[3:24]
            proba_dict[site] = prob
    if "Prob distribution at node " + str(node) + ", by site" in line:
        TAG = 1

file_in.close()

# print(proba_dict)
for site, prob_list in proba_dict.items():
    prob_aa_dict = {}
    if (int(site) > 265) and (int(site) < 408):
        for aa_prob in prob_list:
            prob_aa_dict[aa_prob[0]] = float(aa_prob[2:7])

        sorted_x = sorted(
            prob_aa_dict.items(), key=operator.itemgetter(1), reverse=True
        )

        site_domain = int(site) - 218
        print(
            str(site_domain),
            str(sorted_x[0][0])
            + "("
            + str(round(sorted_x[0][1], 2)).ljust(4, "0")
            + ")",
        )
