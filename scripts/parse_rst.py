#!/usr/bin/python

""" Parse rst """

# Load libraries
import sys

# Load rst file
file_in = open(sys.argv[1], "r")
while 1:
    line = file_in.readline()
    if line == "":
        break
    line = line.rstrip()    # Remove end-of-line
    if line[0:4] == "node":
        tab = line.split("          ")
        print(">"+tab[0].replace(" ", "").replace("#", ""))
        print(tab[1].replace(" ", ""))
file_in.close()
