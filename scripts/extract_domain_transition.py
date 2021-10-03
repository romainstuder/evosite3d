#!/usr/bin/env python3

""" Extract transition between branches """

import sys

# target_branches = ["197..199", "199.200",
#                    "208..209",
#                    "209..329",
#                    "209..210",
#                    "210..265",
#                    "210..211", "211..212", "212..213"]

target_branches = ["208..209", "209..210", "209..265", "208..230"]


# target_position = ["239", "245", "272", "282",
#                    "299", "318", "320", "355",
#                    "366"]

target_position = ['198', '204', '222', '234', '251', '270', '272', '307', '318']


tag = 0
branch = ""
file_in = open(sys.argv[1], "r")
while 1:
    line = file_in.readline()
    if line == "":
        break
    line = line.rstrip()
    if "List of extant and reconstructed sequences" in line:
        tag = 0
    if tag == 1:
        tab = line.split()
        if tab:
            if tab[0] == "Branch":
                branch = tab[2]
            position = tab[0]
            if position in target_position and branch in target_branches:
                print(branch, position, line)
    if "Summary of changes along branches" in line:
        tag = 1

file_in.close()
