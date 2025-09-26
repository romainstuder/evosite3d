import argparse
import sys
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def validate_dataframe(df: pd.DataFrame, required_columns: list) -> None:
    """Validate that dataframe contains required columns."""
    missing_cols = [col for col in required_columns if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    if df.empty:
        raise ValueError("DataFrame is empty")


def load_data(results_file: str, significant_file: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load results and significant phosphosites data."""
    try:
        results_df = pd.read_csv(results_file, sep="\t")
        significant_phosphosites = pd.read_csv(significant_file, sep="\t")

        validate_dataframe(results_df, ["log2_fold_change", "q_value"])
        validate_dataframe(
            significant_phosphosites, ["log2_fold_change", "q_value", "protein", "sequence"]
        )

        return results_df, significant_phosphosites

    except FileNotFoundError as e:
        raise FileNotFoundError(f"Input file not found: {e}")
    except pd.errors.EmptyDataError:
        raise ValueError("One or more input files contain no data")


def create_volcano_plot(
    results_df: pd.DataFrame,
    significant_phosphosites: pd.DataFrame,
    fold_change_threshold: float = 1.0,
    q_value_threshold: float = 0.05,
    output_file: str = "volcano_plot_control_vs_egf5.png",
    figsize: Tuple[int, int] = (10, 8),
    dpi: int = 300,
    show_plot: bool = True,
) -> None:
    """Create volcano plot for phosphoproteome analysis."""

    fig, ax = plt.subplots(figsize=figsize)

    # Plot all points in gray
    ax.scatter(
        results_df["log2_fold_change"],
        -np.log10(results_df["q_value"]),
        alpha=0.6,
        s=20,
        c="gray",
        label="Non-significant",
    )

    # Highlight significant points
    sig_up = significant_phosphosites[
        significant_phosphosites["log2_fold_change"] > fold_change_threshold
    ]
    sig_down = significant_phosphosites[
        significant_phosphosites["log2_fold_change"] < -fold_change_threshold
    ]

    if not sig_up.empty:
        ax.scatter(
            sig_up["log2_fold_change"],
            -np.log10(sig_up["q_value"]),
            alpha=0.8,
            s=30,
            c="red",
            label=f"Upregulated ({len(sig_up)})",
        )

    if not sig_down.empty:
        ax.scatter(
            sig_down["log2_fold_change"],
            -np.log10(sig_down["q_value"]),
            alpha=0.8,
            s=30,
            c="blue",
            label=f"Downregulated ({len(sig_down)})",
        )

    # Add threshold lines
    ax.axhline(
        y=-np.log10(q_value_threshold),
        color="black",
        linestyle="--",
        alpha=0.5,
        label=f"Q-value = {q_value_threshold}",
    )
    ax.axvline(x=fold_change_threshold, color="black", linestyle="--", alpha=0.5)
    ax.axvline(x=-fold_change_threshold, color="black", linestyle="--", alpha=0.5)

    # Customize plot
    ax.set_xlabel("Log2 Fold Change (EGF5 vs Control)")
    ax.set_ylabel("-Log10 Q-value")
    ax.set_title("Volcano Plot: EGF5 vs Control Phosphoproteome")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches="tight")

    if show_plot:
        plt.show()
    else:
        plt.close()

    print(f"Volcano plot saved to {output_file}")


def create_heatmap(
    significant_phosphosites: pd.DataFrame,
    top_n: int = 50,
    output_file: str = "heatmap_top_phosphosites.png",
    figsize: Tuple[int, int] = (12, 10),
    dpi: int = 300,
    show_plot: bool = True,
) -> None:
    """Create heatmap of top changed phosphosites."""

    if len(significant_phosphosites) < top_n:
        print(f"Warning: Only {len(significant_phosphosites)} phosphosites available, using all")
        top_n = len(significant_phosphosites)

    top_changed = significant_phosphosites.nlargest(top_n, "log2_fold_change")

    if top_changed.empty:
        raise ValueError("No phosphosites found for heatmap")

    fig, ax = plt.subplots(figsize=figsize)

    # Create labels
    labels = [f"{row['protein']}_{row['sequence']}" for _, row in top_changed.iterrows()]

    # Create heatmap
    sns.heatmap(
        top_changed[["log2_fold_change"]].T,
        annot=True,
        cmap="RdBu_r",
        center=0,
        yticklabels=["Log2 FC"],
        xticklabels=labels,
        cbar_kws={"label": "Log2 Fold Change"},
        ax=ax,
    )

    ax.set_title(f"Top {top_n} Upregulated Phosphosites (EGF5 vs Control)")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches="tight")

    if show_plot:
        plt.show()
    else:
        plt.close()

    print(f"Heatmap saved to {output_file}")


def main(
    results_file: str,
    significant_file: str,
    fold_change_threshold: float = 1.0,
    q_value_threshold: float = 0.05,
    top_n: int = 50,
    output_dir: str = ".",
    dpi: int = 300,
    show_plots: bool = True,
) -> None:
    """Main function to create MS visualizations."""

    # Load data
    results_df, significant_phosphosites = load_data(results_file, significant_file)

    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Generate plots
    volcano_output = output_path / "volcano_plot_control_vs_egf5.png"
    heatmap_output = output_path / "heatmap_top_phosphosites.png"

    create_volcano_plot(
        results_df,
        significant_phosphosites,
        fold_change_threshold=fold_change_threshold,
        q_value_threshold=q_value_threshold,
        output_file=str(volcano_output),
        dpi=dpi,
        show_plot=show_plots,
    )

    create_heatmap(
        significant_phosphosites,
        top_n=top_n,
        output_file=str(heatmap_output),
        dpi=dpi,
        show_plot=show_plots,
    )

    print("All visualizations completed successfully")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create MS phosphoproteome visualizations")
    parser.add_argument("results_file", help="TSV file with all results data")
    parser.add_argument("significant_file", help="TSV file with significant phosphosites")
    parser.add_argument(
        "--fold-change-threshold",
        type=float,
        default=1.0,
        help="Log2 fold change threshold for significance (default: 1.0)",
    )
    parser.add_argument(
        "--q-value-threshold",
        type=float,
        default=0.05,
        help="Q-value threshold for significance (default: 0.05)",
    )
    parser.add_argument(
        "--top-n", type=int, default=50, help="Number of top phosphosites for heatmap (default: 50)"
    )
    parser.add_argument(
        "--output-dir", default=".", help="Output directory for plots (default: current directory)"
    )
    parser.add_argument("--dpi", type=int, default=300, help="DPI for saved plots (default: 300)")
    parser.add_argument("--no-show", action="store_true", help="Don't display plots interactively")

    args = parser.parse_args()

    try:
        main(
            args.results_file,
            args.significant_file,
            fold_change_threshold=args.fold_change_threshold,
            q_value_threshold=args.q_value_threshold,
            top_n=args.top_n,
            output_dir=args.output_dir,
            dpi=args.dpi,
            show_plots=not args.no_show,
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
