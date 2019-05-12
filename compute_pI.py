#!/usr/bin/python

""" Compute pI """

# Load libraries
import sys
from Bio import SeqIO
from Bio.SeqUtils import ProtParam

# Load file
handle = open(sys.argv[1], "rU")
for record in SeqIO.parse(handle, "fasta"):
    sequence = record.seq.tostring()
    sequence = sequence.replace("-", "")
    analysed_protein = ProtParam.ProteinAnalysis(sequence)

    # Compute some properties
    pI = analysed_protein.isoelectric_point()
    MW = analysed_protein.molecular_weight()
    print(record.id, pI)
handle.close()
