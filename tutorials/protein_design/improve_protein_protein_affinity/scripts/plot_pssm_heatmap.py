import argparse
import logging

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Plot PSSM heatmap from FoldX output")
    parser.add_argument("pdb_id", help="PDB identifier (without .pdb extension)")
    parser.add_argument("workdir", help="Working directory")

    args = parser.parse_args()

    pdb_id = args.pdb_id
    workdir = args.workdir

    # Load compiled data
    df = pd.read_csv(f"{workdir}/{pdb_id}_pssm_output.csv", sep="\t")

    # Define colours for groups
    lut = dict(zip(df["part"].unique(), sns.color_palette()))
    col_colors = df[["position", "part"]].drop_duplicates()["part"].map(lut)
    col_colors.index = df[["position"]].drop_duplicates()["position"]

    #  Pivot data and save it
    df_ddg_table = df.pivot(index="mutant", columns="position", values="ddG")
    sns.clustermap(
        df_ddg_table,
        cmap="inferno_r",
        row_cluster=True,
        col_cluster=False,
        col_colors=col_colors,
        vmin=-1.7,
        vmax=-0,
        figsize=(31, 6),
    )
    plt.savefig(
        f"{workdir}/{pdb_id}_Repair_interface_clustermap.png", bbox_inches="tight", pad_inches=0.0
    )
    logger.info(
        f"Heatmap of {pdb_id} interface saved to {workdir}/{pdb_id}_Repair_interface_clustermap.png"
    )


if __name__ == "__main__":
    main()
