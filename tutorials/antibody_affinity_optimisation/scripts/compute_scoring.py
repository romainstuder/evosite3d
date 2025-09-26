import os
import re
import sys

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy.stats import zscore
from utils import extract_sequence_from_pdb


def find_motifs(sequence, motifs):
    results = {}
    for motif in motifs:
        pattern = motif.replace("X", ".")  # X = any amino acid
        matches = [(m.start(), m.group()) for m in re.finditer(pattern, sequence)]
        if matches:
            results[motif] = matches
    return results


def main():
    pdb_id = sys.argv[1]
    target_chain = sys.argv[2]

    df = pd.read_csv(f"./{pdb_id}_multivariants_output.csv", sep="\t")

    variant_dict = {}

    all_variant_file = open("all_variant_file.fasta", "w")
    pdb_file = f"./{pdb_id}_Repair.pdb"
    seq = extract_sequence_from_pdb(pdb_file, target_chain)
    seq = "".join([res[0] for res in seq])
    all_variant_file.write(f">{pdb_id}_Wild_Type\n")
    all_variant_file.write(seq + "\n")

    variant_dict[f"{pdb_id}_Wild_Type"] = seq
    for variant in df["variant"]:
        # print(variant)
        variant_type = "pairs"
        if variant.count("_") == 2:
            variant_type = "triplets"
        pdb_file = f"./{variant_type}_results/{pdb_id}_Repair_{variant}/{pdb_id}_Repair_1.pdb"
        if os.path.exists(pdb_file):
            seq = extract_sequence_from_pdb(pdb_file, target_chain)
            seq = "".join([res[0] for res in seq])
            all_variant_file.write(f">{variant}\n")
            all_variant_file.write(seq + "\n")
            variant_dict[f"{variant}"] = seq
            # break
    all_variant_file.close()

    # Usage
    properties_list = []
    for name, seq in variant_dict.items():
        X = ProteinAnalysis(seq)
        properties_list.append([name, X.instability_index(), X.isoelectric_point(), X.gravy()])

    df_properties = pd.DataFrame(properties_list)
    df_properties.columns = [
        "variant",
        "instability_index",
        "isoelectric_point",
        "gravy",
    ]

    df_binding_affinity = df[["variant", "total energy"]].copy()
    new_df = df_binding_affinity.merge(df_properties, on="variant", how="outer")

    new_df["total energy"] = new_df["total energy"].fillna(0)
    new_df["total energy_norm"] = -zscore(new_df["total energy"])
    new_df["instability_index_norm"] = -zscore(new_df["instability_index"])
    new_df["isoelectric_point_norm"] = -zscore(new_df["isoelectric_point"])
    new_df["gravy_norm"] = -abs(zscore(new_df["gravy"]))

    # Compute final scoring
    new_df["total_score"] = (
        0.7 * new_df["total energy_norm"]
        + 0.1 * new_df["instability_index_norm"]
        + 0.1 * new_df["isoelectric_point_norm"]
        + 0.1 * new_df["gravy_norm"]
    )
    new_df = new_df.sort_values(["total_score"], ascending=False)
    print("Five best:")
    print(new_df.head(n=5)["variant"])

    for variant in new_df["variant"].head(n=5):
        print(f"load triplets_results/9FWW_Repair_{variant}/{pdb_id}_Repair_1.pdb")
    new_df.head(n=100).to_csv(f"{pdb_id}_variants_ranked.csv", sep="\t", index=False)
    print(f"Wrote to {pdb_id}_variants_ranked.csv")


if __name__ == "__main__":
    main()
