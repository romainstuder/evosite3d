import argparse
import sys
from typing import Optional, Tuple

import pandas as pd


def load_phosphosite_data(input_file: str) -> pd.DataFrame:
    """Load phosphosite data from TSV file."""
    try:
        df = pd.read_csv(input_file, sep="\t")
        if df.empty:
            raise ValueError("Input file is empty")
        return df
    except FileNotFoundError:
        raise FileNotFoundError(f"Input file '{input_file}' not found")
    except pd.errors.EmptyDataError:
        raise ValueError("Input file contains no data")


def filter_phosphorylated_peptides(df: pd.DataFrame) -> pd.DataFrame:
    """Filter dataframe for phosphorylated peptides."""
    if "modifications" not in df.columns:
        raise ValueError("Required column 'modifications' not found in data")

    phospho_df = df[df["modifications"].str.contains("Phospho", na=False)]
    if phospho_df.empty:
        raise ValueError("No phosphorylated peptides found in data")

    return phospho_df


def process_localization_scores(
    phospho_df: pd.DataFrame, score_threshold: float = 0.75
) -> pd.DataFrame:
    """Process localization scores and filter by confidence threshold."""
    if "score" not in phospho_df.columns:
        raise ValueError("Required column 'score' not found in data")

    phospho_df = phospho_df.copy()
    phospho_df["localization_score"] = phospho_df["score"]

    high_confidence = phospho_df[phospho_df["localization_score"] > score_threshold]
    if high_confidence.empty:
        raise ValueError(f"No phosphosites found with score > {score_threshold}")

    return high_confidence


def summarize_by_protein(high_confidence_df: pd.DataFrame) -> pd.DataFrame:
    """Group phosphosites by protein and calculate summary statistics."""
    if "protein_accession" not in high_confidence_df.columns:
        raise ValueError("Required column 'protein_accession' not found in data")

    protein_phospho = (
        high_confidence_df.groupby("protein_accession")
        .agg({"sequence": "count", "localization_score": "mean"})
        .rename(columns={"sequence": "phosphosite_count"})
    )

    return protein_phospho


def save_results(protein_phospho: pd.DataFrame, output_file: Optional[str] = None) -> None:
    """Save results to file if specified."""
    if output_file:
        protein_phospho.to_csv(output_file, sep="\t")
        print(f"Results saved to {output_file}")


def main(
    input_file: str, score_threshold: float = 0.75, output_file: Optional[str] = None
) -> Tuple[int, int]:
    """Main processing function."""
    df = load_phosphosite_data(input_file)
    phospho_df = filter_phosphorylated_peptides(df)
    high_confidence = process_localization_scores(phospho_df, score_threshold)
    protein_phospho = summarize_by_protein(high_confidence)

    save_results(protein_phospho, output_file)

    phosphosite_count = len(high_confidence)
    protein_count = len(protein_phospho)

    print(f"Identified {phosphosite_count} high-confidence phosphosites")
    print(f"Across {protein_count} proteins")

    return phosphosite_count, protein_count


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract and analyze phosphosite information")
    parser.add_argument("input_file", help="Input TSV file with phosphosite data")
    parser.add_argument(
        "--score-threshold",
        type=float,
        default=0.75,
        help="Localization score threshold (default: 0.75)",
    )
    parser.add_argument("--output", "-o", help="Output file for results")

    args = parser.parse_args()

    try:
        main(args.input_file, args.score_threshold, args.output)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
