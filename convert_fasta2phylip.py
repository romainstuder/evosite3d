#!/usr/bin/env python3

"""Convert FASTA to PHYLIP"""

import sys
from Bio import SeqIO

print("Convert FASTA to PHYLIP")

infile = sys.argv[1]
outfile = sys.argv[2]

sequence_list = []  # To keep order of sequence

sequence_dict = {}
for record in SeqIO.parse(open(infile, "rU"), "fasta"):
    tab = record.id.split(" ")
    sequence = str(record.seq).replace(" ", "")
    # print sequence, len(sequence)
    sequence_list.append(tab[0])
    sequence_dict[tab[0]] = sequence
    if "U" in sequence:
        print(tab[0])
        sys.exit()

print("Number of sequences:", len(sequence_dict))


# Test length of the alignment:
alignment_length = 0
for gene in sequence_dict:
    if (alignment_length != 0) and (len(sequence_dict[gene]) != alignment_length):
        print("Error in alignment length, exit on error !!!")
        sys.exit()
    else:
        alignment_length = len(sequence_dict[gene])

number_of_seq = len(sequence_dict)
print("Number of sequences:\t"+str(number_of_seq))
print("Alignment length:\t"+str(alignment_length))
print("Ratio =\t"+str(alignment_length/3))

if alignment_length%3 != 0:
    print("Warning: Hum, your alignment didn't code for nucleotides")

# Length of gene id, can be changed by passing a third argument
name_length = 50
if len(sys.argv) > 3:
    name_length = int(sys.argv[3])


# Write alignment in Phylip format
phyfile = open(outfile, "w")
phyfile.write(str(number_of_seq)+"\t"+str(alignment_length)+"\n")
for gene in sequence_list:
    if len(gene) > name_length:
        gene_name = gene[0:name_length].replace(" ", "")
        if gene_name[-1] == "_":
            gene_name = gene_name[0:-1]
        # elif gene_name[-2] == "_":
        #     gene_name = gene_name[0:-2]
    else:
        gene_name = gene
    phyfile.write(gene_name+"  "+sequence_dict[gene]+"\n")
phyfile.close()
