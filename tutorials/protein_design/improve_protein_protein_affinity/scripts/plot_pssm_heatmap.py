import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def main():
    pdb_id = sys.argv[1]

    # Load compiled data
    df = pd.read_csv(f"{pdb_id}_pssm_output.csv", sep="\t")

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
    plt.savefig(f"{pdb_id}_Repair_interface_clustermap.png", bbox_inches="tight", pad_inches=0.0)


if __name__ == "__main__":
    main()
