#!/usr/bin/env python3

"""Plot dN/dS for sites"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


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

    args = parser.parse_args()

    # Read the data
    df = pd.read_csv(args.input_file, sep=r"\s+", header=None)

    # Create the beb column based on V13 values and threshold
    df["beb"] = df.iloc[:, 12].apply(lambda x: "Yes" if x > args.threshold else "No")

    # Set up the plotting style
    sns.set_style("whitegrid")
    plt.figure(figsize=(args.width, args.height))

    # Create the scatter plot
    ax = sns.scatterplot(
        data=df,
        x=df.columns[0],  # V1 equivalent
        y=df.columns[14],  # V15 equivalent
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
    plt.savefig(args.output, dpi=args.dpi, bbox_inches="tight")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
