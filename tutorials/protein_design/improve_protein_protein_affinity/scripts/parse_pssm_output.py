import argparse
import logging
from itertools import repeat

import numpy as np
import pandas as pd
from utils import AA_LIST, CDR1, CDR2, CDR3, D3TO1, extract_sequence_from_pdb

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Parse PSSM output from FoldX")
    parser.add_argument("pdb_id", help="PDB identifier (without .pdb extension)")
    parser.add_argument("target_chain", help="Target chain identifier")

    args = parser.parse_args()

    pdb_id = args.pdb_id
    target_chain = args.target_chain

    pdb_file = f"{pdb_id}_Repair.pdb"
    seq = extract_sequence_from_pdb(pdb_file, target_chain)

    df_list = []
    for res in seq:
        # print(res)
        # res_to_mutate = "".join(list(res))
        # position = res[2].zfill(3)
        average_file = (
            f"./pssm_results/{pdb_id}_Repair_{res[2].zfill(3)}/Interaction_{pdb_id}_Repair_AC.fxout"
        )
        df = pd.read_csv(average_file, sep="\t")
        df["position"] = int(res[2])
        df["wt"] = res[0]
        df["chain"] = res[1]
        df["mut"] = [x for item in [D3TO1[aa.upper()] for aa in AA_LIST] for x in repeat(item, 2)]
        df["mutant"] = [x for item in AA_LIST for x in repeat(item, 2)]
        df_list.append(df)

    df = pd.concat(df_list)

    df["part"] = "framework"
    df["part"] = np.where(df["position"].isin(CDR1), "CDR1", df["part"])
    df["part"] = np.where(df["position"].isin(CDR2), "CDR2", df["part"])
    df["part"] = np.where(df["position"].isin(CDR3), "CDR3", df["part"])

    df["type"] = np.where(df["Pdb"].str.contains("/WT_"), "wt_dG", "mut_dG")

    df_reduce = df[
        ["position", "wt", "chain", "mut", "mutant", "part", "type", "Interaction Energy"]
    ].copy()
    df_reduce = df_reduce.pivot(
        index=["position", "wt", "chain", "mut", "mutant", "part"],
        columns="type",
        values="Interaction Energy",
    ).reset_index()
    df_reduce["ddG"] = df_reduce["mut_dG"] - df_reduce["wt_dG"]
    df_reduce["mutant_code"] = (
        df_reduce["wt"] + df_reduce["chain"] + df_reduce["position"].astype(str) + df_reduce["mut"]
    )
    df_reduce = df_reduce.sort_values(["ddG"])

    df_reduce.to_csv(f"./{pdb_id}_pssm_output.csv", sep="\t", index=False)
    logger.info(f"output written to ./{pdb_id}_pssm_output.csv")


if __name__ == "__main__":
    main()
