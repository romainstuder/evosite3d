import os
import sys

from utils import extract_sequence_from_pdb


def main():
    pdb_id = sys.argv[1]
    other_chain = sys.argv[2]
    target_chain = sys.argv[3]

    pdb_file = f"{pdb_id}_Repair.pdb"
    seq = extract_sequence_from_pdb(pdb_file, target_chain)

    os.makedirs("./pssm_results/", exist_ok=True)
    command_file = open("command_list.txt", "w")
    for res in seq:
        res_to_mutate = "".join(list(res))
        output_dir = f"./pssm_results/{pdb_id}_Repair_{res[2].zfill(3)}"
        os.makedirs(f"{output_dir}", exist_ok=True)
        command = [
            "./foldx",
            "--command=Pssm",
            f"--analyseComplexChains={other_chain},{target_chain}",
            f"--pdb={pdb_id}_Repair.pdb",
            f"--positions={res_to_mutate}a",
            f"--output-dir={output_dir}",
        ]
        command_file.write(" ".join(command) + "\n")
    command_file.close()

    print("cat command_list.txt | xargs -P 7 -I {} sh -c '{}'")


if __name__ == "__main__":
    main()
