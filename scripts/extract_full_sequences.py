#!/usr/bin/env python3

""" Extract full sequences"""

import sys

tag = 0
branch = ""
file_in = open(sys.argv[1], "r")
while 1:
    line = file_in.readline()
    if line == "":
        break
    line = line.rstrip()
    if "Overall accuracy of the 195 ancestral sequences" in line:
        tag = 0
    if tag == 1 and line:
        if line[0] != " ":
            line = line.replace("e #", "e_")
            tab = line.split()
            if len(tab) > 3:
                # print(tab)
                name = tab[0]
                seq = "".join(tab[1:])
                print(">"+name)
                print(seq)
                # print(tab)
                # site = tab[0]
                # prob = tab[3:24]
                # print(site, prob)
                # print(len(prob))
    if "List of extant and reconstructed sequences" in line:
        tag = 1

file_in.close()
