#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import List, Tuple

import pandas as pd
from analyse_sequences import analyse_sequences
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class SequenceSelector:
    """Class for selecting sequences based on biophysical properties"""

    def __init__(
        self,
        charge_range: Tuple[float, float] = (-2, 5),
        min_hydrophobicity: float = 0.0,
        min_helix_propensity: float = 1.0,
        max_aromatic_content: float = 0.2,
        top_n: int = 10,
    ):
        """Initialize selection criteria"""
        self.charge_range = charge_range
        self.min_hydrophobicity = min_hydrophobicity
        self.min_helix_propensity = min_helix_propensity
        self.max_aromatic_content = max_aromatic_content
        self.top_n = top_n

    def filter_sequences(self, df: pd.DataFrame) -> pd.DataFrame:
        """Filter sequences based on defined criteria"""
        filtered = df[
            (df["charge"].between(self.charge_range[0], self.charge_range[1]))
            & (df["hydrophobicity"] > self.min_hydrophobicity)
            & (df["helix_propensity"] > self.min_helix_propensity)
            & (df["aromatic_content"] < self.max_aromatic_content)
        ]

        filtered = filtered.sort_values("helix_propensity", ascending=False)
        return filtered.head(self.top_n)

    def write_filtered_sequences(
        self, top_sequences: pd.DataFrame, fasta_file: Path, output_file: Path
    ) -> List[SeqRecord]:
        """Write filtered sequences to output file"""
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        filtered_seqs = [seq for seq in sequences if seq.id in top_sequences["id"].values]

        SeqIO.write(filtered_seqs, output_file, "fasta")
        return filtered_seqs

    def select_sequences(self, fasta_file: Path, output_file: Path = None) -> pd.DataFrame:
        """Main method to analyze and select sequences"""
        if output_file is None:
            output_file = fasta_file.parent / "filtered_sequences.fasta"

        df = analyse_sequences(str(fasta_file))
        top_sequences = self.filter_sequences(df)

        if not top_sequences.empty:
            self.write_filtered_sequences(top_sequences, fasta_file, output_file)
            print(f"Selected {len(top_sequences)} sequences, written to {output_file}")
        else:
            print("No sequences met the filtering criteria")

        return top_sequences


def main():
    """Main function with command line interface"""
    parser = argparse.ArgumentParser(description="Select sequences based on biophysical properties")
    parser.add_argument("--input_file", type=str, help="Input FASTA file")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="filtered_sequences.fasta",
        help="Output FASTA file (default: filtered_sequences.fasta)",
    )
    parser.add_argument("--charge-min", type=float, default=-2, help="Minimum charge (default: -2)")
    parser.add_argument("--charge-max", type=float, default=5, help="Maximum charge (default: 5)")
    parser.add_argument(
        "--min-hydrophobicity",
        type=float,
        default=0.0,
        help="Minimum hydrophobicity (default: 0.0)",
    )
    parser.add_argument(
        "--min-helix-propensity",
        type=float,
        default=1.0,
        help="Minimum helix propensity (default: 1.0)",
    )
    parser.add_argument(
        "--max-aromatic", type=float, default=0.2, help="Maximum aromatic content (default: 0.2)"
    )
    parser.add_argument(
        "-n",
        "--top-n",
        type=int,
        default=10,
        help="Number of top sequences to select (default: 10)",
    )
    parser.add_argument("--verbose", action="store_true", help="Print detailed results")

    args = parser.parse_args()

    input_path = Path(args.input_file)
    output_path = Path(args.output)

    if not input_path.exists():
        print(f"Error: Input file {input_path} does not exist")
        return

    selector = SequenceSelector(
        charge_range=(args.charge_min, args.charge_max),
        min_hydrophobicity=args.min_hydrophobicity,
        min_helix_propensity=args.min_helix_propensity,
        max_aromatic_content=args.max_aromatic,
        top_n=args.top_n,
    )

    top_seqs = selector.select_sequences(input_path, output_path)

    if args.verbose and not top_seqs.empty:
        print("\nSelected sequences:")
        print(top_seqs[["id", "charge", "hydrophobicity", "helix_propensity", "aromatic_content"]])


if __name__ == "__main__":
    main()
