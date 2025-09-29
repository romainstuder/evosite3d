import os
import sys
from itertools import combinations

import pandas as pd


def unique_pairs_only(lst):
    """Get unique pairs without self-pairing"""
    return list(combinations(lst, 2))


def unique_triplets_only(lst):
    """Get unique pairs without self-pairing"""
    return list(combinations(lst, 3))


def remove_same_positions(lst):
    new_list = []
    for pair in lst:
        pos_list = [pos[2:-1] for pos in pair]
        if len(set(pos_list)) == len(pos_list):
            new_list.append(pair)
    return new_list


def write_mutant_pairs(df_best, pdb_id):
    # Build commands for pairwise
    output_pairs_results_folder = "pairs_results"
    os.makedirs(output_pairs_results_folder, exist_ok=True)

    command_file = open("command_list_pairs.txt", "w")

    unique_pairs_list = unique_pairs_only(df_best["mutant_code"])
    unique_pairs_list = remove_same_positions(unique_pairs_list)
    print(f"Number of variants (pairs): {len(unique_pairs_list)}")

    for pair in unique_pairs_list:
        # print(pair)
        res_to_mutate = ",".join(list(pair))
        # print(res_to_mutate)
        # position = res[2].zfill(3)
        output_dir = (
            f"{output_pairs_results_folder}/{pdb_id}_Repair_{res_to_mutate.replace(',', '_')}"
        )
        os.makedirs(f"{output_dir}", exist_ok=True)
        with open(f"{output_dir}/individual_list.txt", "w") as f:
            f.write(f"{res_to_mutate};\n")

        #  Write command
        command = [
            "./foldx",
            "--command=BuildModel",
            f"--pdb={pdb_id}_Repair.pdb",
            f"--mutant-file={output_dir}/individual_list.txt",
            f"--output-dir={output_dir}",
        ]
        #  print(" ".join(command))
        command_file.write(" ".join(command) + "\n")
    command_file.close()

    print("cat command_list_pairs.txt | xargs -P 7 -I {} sh -c '{}'")


def write_mutant_triplets(df_best, pdb_id):
    # Build commands for triplets
    unique_triplets_list = unique_triplets_only(
        df_best.sort_values(["position", "wt"])["mutant_code"]
    )
    unique_triplets_list = remove_same_positions(unique_triplets_list)
    print(f"Number of variants (triplets): {len(unique_triplets_list)}")
    output_triplets_results_folder = "triplets_results"

    os.makedirs(output_triplets_results_folder, exist_ok=True)
    command_file = open("command_list_triplets.txt", "w")
    for pair in unique_triplets_list:
        # print(pair)
        res_to_mutate = ",".join(list(pair))
        # print(res_to_mutate)
        # position = res[2].zfill(3)
        output_dir = (
            f"{output_triplets_results_folder}/{pdb_id}_Repair_{res_to_mutate.replace(',', '_')}"
        )
        os.makedirs(output_dir, exist_ok=True)
        with open(f"{output_dir}/individual_list.txt", "w") as f:
            f.write(f"{res_to_mutate};\n")

        # Write command
        command = [
            "./foldx",
            "--command=BuildModel",
            f"--pdb={pdb_id}_Repair.pdb",
            f"--mutant-file={output_dir}/individual_list.txt",
            f"--output-dir={output_dir}",
        ]
        # print(" ".join(command))
        command_file.write(" ".join(command) + "\n")
    command_file.close()

    print("cat command_list_triplets.txt | xargs -P 7 -I {} sh -c '{}'")


def main():
    pdb_id = sys.argv[1]
    dgg_threshold = float(sys.argv[2])

    df = pd.read_csv(f"{pdb_id}_pssm_output.csv", sep="\t")
    df_best = df[df["ddG"] <= dgg_threshold].copy()

    print(df_best.head())

    write_mutant_pairs(df_best, pdb_id)
    write_mutant_triplets(df_best, pdb_id)


if __name__ == "__main__":
    main()
