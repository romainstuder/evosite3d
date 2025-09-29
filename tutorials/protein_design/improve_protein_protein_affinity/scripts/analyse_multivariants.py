import argparse
import logging

import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Analyse multivariants from FoldX output")
    parser.add_argument("pdb_id", help="PDB identifier (without .pdb extension)")

    args = parser.parse_args()

    pdb_id = args.pdb_id

    df_list = []
    for variants_type in ["pairs", "triplets"]:
        # Load variants_type
        file = open(f"./command_list_{variants_type}.txt", "r")
        line = file.readline()
        folder_list = []
        while line:
            line = file.readline().rstrip()
            if line != "":
                output_folder = line.split("--output-dir=")[-1]
                folder_list.append(output_folder)
        file.close()
        logger.info(f"Number of {variants_type} variants: {len(folder_list)}")

        for output_folder in folder_list:
            average_file = f"{output_folder}/Average_{pdb_id}_Repair.fxout"
            df = pd.read_csv(average_file, sep="\t", skiprows=8)
            # print(df)
            df["variant"] = output_folder.replace(f"{variants_type}_results/{pdb_id}_Repair_", "")
            # print(output_folder.replace(f"{variants_type}_results/{pdb_id}_Repair_", ""))
            df_list.append(df)

    df_output = pd.concat(df_list)
    # assert len(df_output) == len(df_list), "Number of output files and number of input files is
    # different"
    df_output = df_output.sort_values(["total energy"])
    cols = ["variant"] + [col for col in df_output.columns if col != "variant"]
    df_output = df_output.reindex(columns=cols)
    df_output.to_csv(f"./{pdb_id}_multivariants_output.csv", index=False, sep="\t")


if __name__ == "__main__":
    main()
