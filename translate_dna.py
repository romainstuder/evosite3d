#!/usr/bin/env python3

"""Translate DNA"""


import sys
from Bio import SeqIO

handle = open(sys.argv[1], "rU")
for record in SeqIO.parse(handle, "fasta"):
    print(">"+record.id)
    print(record.seq.translate(to_stop=True))
handle.close()
