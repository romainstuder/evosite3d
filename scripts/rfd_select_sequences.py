#!/usr/bin/env python3

from Bio import SeqIO
from rfd_analyse_sequences import analyse_sequences


def filter_sequences(df, fasta_file):
    """Filter sequences based on desired properties"""
    # Define criteria for two-helix bundle
    filtered = df[
        (df["charge"].between(-2, 5))
        & (df["hydrophobicity"] > 0)  # Slightly positive
        & (df["helix_propensity"] > 1.0)  # Hydrophobic core
        & (df["aromatic_content"] < 0.2)  # High helix propensity  # Not too many aromatics
    ]

    # Sort by helix propensity
    filtered = filtered.sort_values("helix_propensity", ascending=False)

    # Get top sequences
    top_sequences = filtered.head(10)

    # Write filtered sequences
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    filtered_seqs = []

    for seq in sequences:
        if seq.id in top_sequences["id"].values:
            filtered_seqs.append(seq)

    SeqIO.write(filtered_seqs, "filtered_sequences.fasta", "fasta")

    return top_sequences


# First, analyze the sequences to create the DataFrame
df = analyse_sequences("./output/sequences/sequences.fasta")

# Then filter and select
top_seqs = filter_sequences(df, "./output/sequences/sequences.fasta")
print("Top 10 sequences selected")
print(top_seqs[["id", "charge", "hydrophobicity", "helix_propensity"]])
