#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO


def analyse_sequences(fasta_file):
    """Analyse properties of designed sequences"""
    sequences = list(SeqIO.parse(fasta_file, "fasta"))

    results = []
    for seq_record in sequences:
        seq = str(seq_record.seq)

        # Calculate properties
        properties = {
            "id": seq_record.id,
            "length": len(seq),
            "charge": calculate_charge(seq),
            "hydrophobicity": calculate_hydrophobicity(seq),
            "helix_propensity": calculate_helix_propensity(seq),
            "aromatic_content": sum(seq.count(aa) for aa in "FWY") / len(seq),
        }

        results.append(properties)

    return pd.DataFrame(results)


def calculate_charge(sequence, pH=7.0):
    """Calculate net charge at given pH"""
    positive = sequence.count("K") + sequence.count("R") + sequence.count("H") * 0.1
    negative = sequence.count("D") + sequence.count("E")
    return positive - negative


def calculate_hydrophobicity(sequence):
    """Calculate average hydrophobicity (Kyte-Doolittle scale)"""
    kd_scale = {
        "I": 4.5,
        "V": 4.2,
        "L": 3.8,
        "F": 2.8,
        "C": 2.5,
        "M": 1.9,
        "A": 1.8,
        "G": -0.4,
        "T": -0.7,
        "S": -0.8,
        "W": -0.9,
        "Y": -1.3,
        "P": -1.6,
        "H": -3.2,
        "E": -3.5,
        "Q": -3.5,
        "D": -3.5,
        "N": -3.5,
        "K": -3.9,
        "R": -4.5,
    }
    return sum(kd_scale.get(aa, 0) for aa in sequence) / len(sequence)


def calculate_helix_propensity(sequence):
    """Calculate average helix propensity"""
    helix_scale = {
        "A": 1.42,
        "L": 1.21,
        "M": 1.20,
        "E": 1.17,
        "K": 1.16,
        "F": 1.13,
        "Q": 1.10,
        "I": 1.08,
        "W": 1.08,
        "V": 1.06,
        "D": 1.01,
        "H": 1.00,
        "R": 0.98,
        "T": 0.87,
        "S": 0.79,
        "C": 0.79,
        "Y": 0.69,
        "N": 0.67,
        "P": 0.57,
        "G": 0.57,
    }
    return sum(helix_scale.get(aa, 1.0) for aa in sequence) / len(sequence)


def main():
    """Main function to analyse sequences"""
    import argparse

    parser = argparse.ArgumentParser(description="Analyse properties of designed sequences")
    parser.add_argument("--fasta_file", help="Input FASTA file containing sequences")
    parser.add_argument(
        "-o",
        "--output",
        default="sequence_properties.png",
        help="Output plot file (default: sequence_properties.png)",
    )

    args = parser.parse_args()

    # Analyse designed sequences
    df = analyse_sequences(args.fasta_file)
    print(df.describe())

    # Plot distributions
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    df["charge"].hist(ax=axes[0, 0], bins=20)
    axes[0, 0].set_title("Charge Distribution")
    df["hydrophobicity"].hist(ax=axes[0, 1], bins=20)
    axes[0, 1].set_title("Hydrophobicity Distribution")
    df["helix_propensity"].hist(ax=axes[1, 0], bins=20)
    axes[1, 0].set_title("Helix Propensity Distribution")
    df["aromatic_content"].hist(ax=axes[1, 1], bins=20)
    axes[1, 1].set_title("Aromatic Content Distribution")
    plt.tight_layout()
    plt.savefig(args.output)


if __name__ == "__main__":
    main()
