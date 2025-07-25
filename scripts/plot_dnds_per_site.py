#!/usr/bin/env python3

"""Plot dN/dS for sites with selective pressure analysis"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def load_and_validate_data(
    input_file: str, position_col: int = 0, beb_col: int = 12, dnds_col: int = 14
) -> pd.DataFrame:
    """Load and validate input data file."""
    try:
        input_path = Path(input_file)
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")

        if not input_path.is_file():
            raise ValueError(f"Path is not a file: {input_file}")

        # Read the data
        df = pd.read_csv(input_file, sep=r"\s+", header=None)

        if df.empty:
            raise ValueError("Input file is empty")

        # Validate required columns exist
        max_col = max(position_col, beb_col, dnds_col)
        if len(df.columns) <= max_col:
            raise ValueError(
                f"File has only {len(df.columns)} columns, but column {max_col + 1} is required"
            )

        # Check for numeric data in required columns
        numeric_cols = [position_col, beb_col, dnds_col]
        for col in numeric_cols:
            if not pd.api.types.is_numeric_dtype(df.iloc[:, col]):
                try:
                    df.iloc[:, col] = pd.to_numeric(df.iloc[:, col], errors="coerce")
                    if df.iloc[:, col].isna().any():
                        raise ValueError(f"Column {col + 1} contains non-numeric values")
                except Exception:
                    raise ValueError(f"Column {col + 1} cannot be converted to numeric")

        return df

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error reading file: {e}", file=sys.stderr)
        sys.exit(1)


def create_plot(
    df: pd.DataFrame,
    threshold: float,
    position_col: int,
    beb_col: int,
    dnds_col: int,
    width: float,
    height: float,
    output_file: str,
    dpi: int,
    show_plot: bool,
) -> None:
    """Create and save the dN/dS scatter plot."""
    try:
        # Create the beb column based on threshold
        df["beb"] = df.iloc[:, beb_col].apply(lambda x: "Yes" if x > threshold else "No")

        # Set up the plotting style
        sns.set_style("whitegrid")
        plt.figure(figsize=(width, height))

        # Create the scatter plot
        ax = sns.scatterplot(
            data=df,
            x=df.columns[position_col],
            y=df.columns[dnds_col],
            hue="beb",
            palette={"No": "black", "Yes": "red"},
            legend=False,
        )

        # Add horizontal line at y=1
        plt.axhline(y=1, color="black", linestyle="-", linewidth=0.8)

        # Customize the plot
        plt.xlabel("Residue position")
        plt.ylabel("Selective pressure [dN/dS]")

        # Remove grid and customize background
        sns.despine()
        ax.grid(False)

        # Save the plot
        plt.tight_layout()
        plt.savefig(output_file, dpi=dpi, bbox_inches="tight")

        print(f"Plot saved to: {output_file}")

        if show_plot:
            plt.show()
        else:
            plt.close()

    except Exception as e:
        print(f"Error creating plot: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Generate BEB scatter plot with selective pressure analysis"
    )
    parser.add_argument("input_file", help="Input data file (whitespace-separated)")
    parser.add_argument(
        "-o", "--output", default="beb.png", help="Output filename (default: beb.png)"
    )
    parser.add_argument(
        "-t", "--threshold", type=float, default=0.50, help="BEB threshold value (default: 0.50)"
    )
    parser.add_argument(
        "--width", type=float, default=4, help="Figure width in inches (default: 4)"
    )
    parser.add_argument(
        "--height", type=float, default=3, help="Figure height in inches (default: 3)"
    )
    parser.add_argument("--dpi", type=int, default=300, help="DPI for saved figure (default: 300)")
    parser.add_argument("--no-show", action="store_true", help="Don't display the plot")
    parser.add_argument(
        "--position-col",
        type=int,
        default=1,
        help="Column number for position data (1-based, default: 1)",
    )
    parser.add_argument(
        "--beb-col",
        type=int,
        default=13,
        help="Column number for BEB values (1-based, default: 13)",
    )
    parser.add_argument(
        "--dnds-col",
        type=int,
        default=15,
        help="Column number for dN/dS values (1-based, default: 15)",
    )

    args = parser.parse_args()

    # Validate arguments
    if args.threshold < 0 or args.threshold > 1:
        print("Error: Threshold must be between 0 and 1", file=sys.stderr)
        sys.exit(1)

    if args.width <= 0 or args.height <= 0:
        print("Error: Width and height must be positive", file=sys.stderr)
        sys.exit(1)

    if args.dpi <= 0:
        print("Error: DPI must be positive", file=sys.stderr)
        sys.exit(1)

    # Convert 1-based column numbers to 0-based indices
    position_col = args.position_col - 1
    beb_col = args.beb_col - 1
    dnds_col = args.dnds_col - 1

    if min(position_col, beb_col, dnds_col) < 0:
        print("Error: Column numbers must be positive", file=sys.stderr)
        sys.exit(1)

    # Load and validate data
    df = load_and_validate_data(args.input_file, position_col, beb_col, dnds_col)

    # Create and save plot
    create_plot(
        df,
        args.threshold,
        position_col,
        beb_col,
        dnds_col,
        args.width,
        args.height,
        args.output,
        args.dpi,
        not args.no_show,
    )


if __name__ == "__main__":
    main()
