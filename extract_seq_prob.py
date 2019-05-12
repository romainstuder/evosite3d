#!/usr/bin/env python3

""" Extract probabilityes """

import sys

tag = 0
branch = ""
file_in = open(sys.argv[1], "r")
while 1:
    line = file_in.readline()
    if line == "":
        break
    line = line.rstrip()
    if "Prob distribution at node 210, by site" in line:
        tag = 0
    if tag == 1:
        tab = line.split()
        if len(tab) > 3:
            # print(tab)
            site = tab[0]
            prob = tab[3:24]
            print(site, prob)
            print(len(prob))
    if "Prob distribution at node 209, by site" in line:
        tag = 1

file_in.close()
