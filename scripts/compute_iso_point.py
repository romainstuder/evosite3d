#!/usr/bin/env python3

"""Compute the isoelectric point (pI) and molecular weight (MW) of proteins from a FASTA file."""

import argparse
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqUtils import ProtParam


def analyze_protein_sequence(sequence):
    clean_sequence = sequence.replace("-", "")
    if not clean_sequence:
        raise ValueError("Empty sequence after cleaning")

    try:
        protein = ProtParam.ProteinAnalysis(clean_sequence)
        return protein.isoelectric_point(), protein.molecular_weight()
    except Exception as e:
        raise ValueError(f"Failed to analyze protein sequence: {e}")


def validate_fasta_file(fasta_file):
    file_path = Path(fasta_file)
    if not file_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
    if not file_path.is_file():
        raise ValueError(f"Path is not a file: {fasta_file}")


def compute_properties(fasta_file):
    try:
        validate_fasta_file(fasta_file)

        print("Sequence_ID\tpI\tMW(Da)")

        with open(fasta_file, "r") as handle:
            records_processed = 0
            for record in SeqIO.parse(handle, "fasta"):
                try:
                    pI, MW = analyze_protein_sequence(str(record.seq))
                    print(f"{record.id}\t{pI:.2f}\t{MW:.2f}")
                    records_processed += 1
                except ValueError as e:
                    print(f"Warning: Skipping {record.id}: {e}", file=sys.stderr)
                    continue

            if records_processed == 0:
                raise ValueError("No valid protein sequences found in the file")

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Compute pI and MW from a FASTA file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Output format: Sequence_ID\\tpI\\tMW(Da)",
    )
    parser.add_argument("fasta_file", help="Input FASTA file with protein sequences.")
    args = parser.parse_args()

    compute_properties(args.fasta_file)


if __name__ == "__main__":
    main()
