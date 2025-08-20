#!/usr/bin/env python3
import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def validate_files(files: List[str]) -> None:
    """Validate that all input files exist."""
    missing_files = [f for f in files if not Path(f).exists()]
    if missing_files:
        raise FileNotFoundError(f"Missing files: {missing_files}")


def load_and_process_data(files: List[str], condition: str) -> pd.DataFrame:
    """Load and combine data from multiple files with error handling."""
    validate_files(files)

    dfs = []
    for file_path in files:
        try:
            df = pd.read_csv(file_path, sep="\t")

            # Validate required columns
            required_cols = [
                "modifications",
                "localization_score",
                "intensity",
                "sequence",
                "protein_accession",
            ]
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                raise ValueError(f"Missing columns in {file_path}: {missing_cols}")

            # Extract replicate number from filename
            try:
                replicate = file_path.split("_rep")[1].split("_")[0]
            except IndexError:
                raise ValueError(f"Cannot extract replicate number from filename: {file_path}")

            df["condition"] = condition
            df["replicate"] = replicate
            df["source_file"] = file_path
            dfs.append(df)

        except Exception as e:
            raise ValueError(f"Error processing file {file_path}: {e}")

    if not dfs:
        raise ValueError(f"No valid data loaded for condition {condition}")

    return pd.concat(dfs, ignore_index=True)


def filter_phosphosites(
    data: pd.DataFrame, localization_threshold: float = 0.75, min_intensity: float = 0.0
) -> pd.DataFrame:
    """Filter for high-confidence phosphosites."""
    if data.empty:
        raise ValueError("Input data is empty")

    # Filter for phosphorylated peptides
    phospho_mask = data["modifications"].str.contains("Phospho", na=False)
    if not phospho_mask.any():
        raise ValueError("No phosphorylated peptides found in data")

    # Apply filters
    filtered_data = data[
        phospho_mask
        & (data["localization_score"] > localization_threshold)
        & (data["intensity"] > min_intensity)
    ]

    if filtered_data.empty:
        raise ValueError(
            f"No phosphosites pass filters (localization > {localization_threshold}, "
            f"intensity > {min_intensity})"
        )

    return filtered_data


def create_intensity_matrix(phospho_data: pd.DataFrame) -> pd.DataFrame:
    """Create intensity matrix for statistical analysis."""
    intensity_matrix = phospho_data.pivot_table(
        index=["sequence", "protein_accession", "modifications"],
        columns=["condition", "replicate"],
        values="intensity",
        fill_value=0,
    )

    if intensity_matrix.empty:
        raise ValueError("No data available for creating intensity matrix")

    return intensity_matrix


def perform_statistical_analysis(
    intensity_matrix: pd.DataFrame,
    condition1: str = "Control",
    condition2: str = "EGF5",
    min_replicates: int = 2,
    statistical_test: str = "ttest",
) -> pd.DataFrame:
    """Perform statistical testing on phosphosite data."""
    results = []

    for phosphosite in intensity_matrix.index:
        try:
            # Get intensities for each condition
            cond1_intensities = intensity_matrix.loc[phosphosite, condition1].values
            cond2_intensities = intensity_matrix.loc[phosphosite, condition2].values

            # Remove zeros for statistical testing
            cond1_nonzero = cond1_intensities[cond1_intensities > 0]
            cond2_nonzero = cond2_intensities[cond2_intensities > 0]

            # Check minimum replicate requirement
            if len(cond1_nonzero) >= min_replicates and len(cond2_nonzero) >= min_replicates:
                # Log2 transformation
                log2_cond1 = np.log2(cond1_nonzero)
                log2_cond2 = np.log2(cond2_nonzero)

                # Statistical test
                if statistical_test == "ttest":
                    t_stat, p_value = stats.ttest_ind(log2_cond1, log2_cond2)
                elif statistical_test == "mannwhitney":
                    t_stat, p_value = stats.mannwhitneyu(
                        log2_cond1, log2_cond2, alternative="two-sided"
                    )
                else:
                    raise ValueError(f"Unknown statistical test: {statistical_test}")

                # Calculate fold change
                mean_cond1 = np.mean(log2_cond1)
                mean_cond2 = np.mean(log2_cond2)
                log2_fc = mean_cond2 - mean_cond1

                results.append(
                    {
                        "sequence": phosphosite[0],
                        "protein": phosphosite[1],
                        "modification": phosphosite[2],
                        "log2_fold_change": log2_fc,
                        "p_value": p_value,
                        "test_statistic": t_stat,
                        f"mean_{condition1.lower()}": mean_cond1,
                        f"mean_{condition2.lower()}": mean_cond2,
                        f"n_{condition1.lower()}": len(cond1_nonzero),
                        f"n_{condition2.lower()}": len(cond2_nonzero),
                    }
                )
        except Exception as e:
            print(f"Warning: Error analyzing phosphosite {phosphosite}: {e}")
            continue

    if not results:
        raise ValueError("No valid statistical results generated")

    return pd.DataFrame(results)


def apply_multiple_testing_correction(
    results_df: pd.DataFrame, method: str = "fdr_bh"
) -> pd.DataFrame:
    """Apply multiple testing correction."""
    if "p_value" not in results_df.columns:
        raise ValueError("p_value column not found in results")

    if results_df["p_value"].isna().any():
        raise ValueError("NaN values found in p_value column")

    _, q_values, _, _ = multipletests(results_df["p_value"], method=method)
    results_df["q_value"] = q_values
    results_df["correction_method"] = method

    return results_df


def filter_significant_results(
    results_df: pd.DataFrame, q_threshold: float = 0.05, fc_threshold: float = 1.0
) -> pd.DataFrame:
    """Filter for significant results."""
    significant = results_df[
        (results_df["q_value"] < q_threshold)
        & (np.abs(results_df["log2_fold_change"]) > fc_threshold)
    ]

    return significant


def generate_summary_statistics(
    results_df: pd.DataFrame, significant_df: pd.DataFrame, fc_threshold: float = 1.0
) -> Dict[str, Any]:
    """Generate summary statistics."""
    upregulated = significant_df[significant_df["log2_fold_change"] > fc_threshold]
    downregulated = significant_df[significant_df["log2_fold_change"] < -fc_threshold]

    summary = {
        "total_phosphosites_analyzed": len(results_df),
        "significant_phosphosites": len(significant_df),
        "upregulated": len(upregulated),
        "downregulated": len(downregulated),
        "fold_change_threshold": fc_threshold,
        "q_value_threshold": significant_df["q_value"].max() if not significant_df.empty else None,
    }

    return summary


def save_results(
    results_df: pd.DataFrame,
    significant_df: pd.DataFrame,
    summary_stats: Dict[str, Any],
    output_prefix: str = "comparative_analysis",
) -> None:
    """Save analysis results to files."""
    # Save all results
    all_results_file = f"{output_prefix}_all_results.tsv"
    results_df.to_csv(all_results_file, sep="\t", index=False)
    print(f"All results saved to: {all_results_file}")

    # Save significant results
    sig_results_file = f"{output_prefix}_significant_phosphosites.tsv"
    significant_df.to_csv(sig_results_file, sep="\t", index=False)
    print(f"Significant results saved to: {sig_results_file}")

    # Save summary statistics
    summary_file = f"{output_prefix}_summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary_stats, f, indent=2)
    print(f"Summary statistics saved to: {summary_file}")


def main(
    control_files: List[str],
    treatment_files: List[str],
    control_name: str = "Control",
    treatment_name: str = "Treatment",
    localization_threshold: float = 0.75,
    min_intensity: float = 0.0,
    q_threshold: float = 0.05,
    fc_threshold: float = 1.0,
    min_replicates: int = 2,
    statistical_test: str = "ttest",
    correction_method: str = "fdr_bh",
    output_prefix: str = "comparative_analysis",
) -> None:
    """Main analysis function."""

    print(f"Loading {control_name} data from {len(control_files)} files...")
    control_data = load_and_process_data(control_files, control_name)

    print(f"Loading {treatment_name} data from {len(treatment_files)} files...")
    treatment_data = load_and_process_data(treatment_files, treatment_name)

    # Combine datasets
    print("Combining datasets...")
    all_data = pd.concat([control_data, treatment_data], ignore_index=True)

    # Filter for high-confidence phosphosites
    print(
        f"Filtering phosphosites (localization > {localization_threshold}, "
        f"intensity > {min_intensity})..."
    )
    phospho_data = filter_phosphosites(all_data, localization_threshold, min_intensity)
    print(f"Found {len(phospho_data)} high-confidence phosphosites")

    # Create intensity matrix
    print("Creating intensity matrix...")
    intensity_matrix = create_intensity_matrix(phospho_data)
    print(f"Matrix dimensions: {intensity_matrix.shape}")

    # Perform statistical analysis
    print(f"Performing statistical analysis ({statistical_test})...")
    results_df = perform_statistical_analysis(
        intensity_matrix, control_name, treatment_name, min_replicates, statistical_test
    )
    print(f"Analyzed {len(results_df)} phosphosites")

    # Apply multiple testing correction
    print(f"Applying multiple testing correction ({correction_method})...")
    results_df = apply_multiple_testing_correction(results_df, correction_method)

    # Filter significant results
    print(f"Filtering significant results (q < {q_threshold}, |log2FC| > {fc_threshold})...")
    significant_df = filter_significant_results(results_df, q_threshold, fc_threshold)

    # Generate summary statistics
    summary_stats = generate_summary_statistics(results_df, significant_df, fc_threshold)

    # Print summary
    print("\n=== ANALYSIS SUMMARY ===")
    for key, value in summary_stats.items():
        print(f"{key}: {value}")

    # Save results
    save_results(results_df, significant_df, summary_stats, output_prefix)

    print("\nAnalysis completed successfully!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MS/MS comparative phosphoproteome analysis")
    parser.add_argument(
        "--control-files", nargs="+", required=True, help="List of control condition TSV files"
    )
    parser.add_argument(
        "--treatment-files", nargs="+", required=True, help="List of treatment condition TSV files"
    )
    parser.add_argument(
        "--control-name", default="Control", help="Name for control condition (default: Control)"
    )
    parser.add_argument(
        "--treatment-name",
        default="Treatment",
        help="Name for treatment condition (default: Treatment)",
    )
    parser.add_argument(
        "--localization-threshold",
        type=float,
        default=0.75,
        help="Localization score threshold (default: 0.75)",
    )
    parser.add_argument(
        "--min-intensity",
        type=float,
        default=0.0,
        help="Minimum intensity threshold (default: 0.0)",
    )
    parser.add_argument(
        "--q-threshold",
        type=float,
        default=0.05,
        help="Q-value significance threshold (default: 0.05)",
    )
    parser.add_argument(
        "--fc-threshold", type=float, default=1.0, help="Log2 fold change threshold (default: 1.0)"
    )
    parser.add_argument(
        "--min-replicates",
        type=int,
        default=2,
        help="Minimum replicates per condition (default: 2)",
    )
    parser.add_argument(
        "--statistical-test",
        choices=["ttest", "mannwhitney"],
        default="ttest",
        help="Statistical test to use (default: ttest)",
    )
    parser.add_argument(
        "--correction-method",
        default="fdr_bh",
        help="Multiple testing correction method (default: fdr_bh)",
    )
    parser.add_argument(
        "--output-prefix",
        default="comparative_analysis",
        help="Output files prefix (default: comparative_analysis)",
    )

    args = parser.parse_args()

    try:
        main(
            control_files=args.control_files,
            treatment_files=args.treatment_files,
            control_name=args.control_name,
            treatment_name=args.treatment_name,
            localization_threshold=args.localization_threshold,
            min_intensity=args.min_intensity,
            q_threshold=args.q_threshold,
            fc_threshold=args.fc_threshold,
            min_replicates=args.min_replicates,
            statistical_test=args.statistical_test,
            correction_method=args.correction_method,
            output_prefix=args.output_prefix,
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
