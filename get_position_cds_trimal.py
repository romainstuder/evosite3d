#!/usr/bin/env python2

# Load libraries
import sys
from Bio import SeqIO

# Load untrimmed sequences from full alignment
sequence_dict = {}
handle = open(sys.argv[1], "rU")
for record in SeqIO.parse(handle, "fasta"):
    sequence = record.seq.translate(to_stop=True, gap="-")
    sequence_dict[record.id] = sequence
handle.close()

# Load columns kept after Trimal
column_list = []
file_in = open(sys.argv[2], "r")
while 1:
    line = file_in.readline()
    if line == "":
        break
    line = line.rstrip()
    if line[0] == "#":
        line = line.replace(",", "")
        tab = line.split()
        column_list = tab[1:]
file_in.close()

# Get corresponding positions per column
corresponding_dict = {}
for gene in sequence_dict:
    j = 0
    corr_dict = {}
    for i in range(len(sequence_dict[gene])):
        if sequence_dict[gene][i] != "-":
            j = j+1
            corr_dict[i+1] = j
    corresponding_dict[gene] = corr_dict

# Define reference gene
gene = sys.argv[3]

# Define position from CodeML
site = sys.argv[4]

# Assign column to the full alignment
column_in_full = int(column_list[int(site)*3])/3

# Print output
print site, corresponding_dict[gene][int(column_in_full)], sequence_dict[gene][column_in_full-1]
