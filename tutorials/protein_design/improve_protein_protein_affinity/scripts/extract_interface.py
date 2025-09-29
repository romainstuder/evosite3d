import sys

import pandas as pd


def extract_interface(pdb_id, chain_id):
    interface_file = f"./Interface_Residues_{pdb_id}_Repair_AC.fxout"
    df = pd.read_csv(interface_file, sep="\t", skiprows=10, header=0)
    interface_residues = df.columns.tolist()
    interface_residues.remove(df.columns[-1])

    pos_in_chain = {}

    for res in interface_residues:
        chain = res[1]
        pos = res[2:]
        pos_in_chain.setdefault(chain, []).append(pos)

    return pos_in_chain[chain_id]


def main():
    pdb_id = sys.argv[1]
    chain_id = sys.argv[2]

    residue_list = extract_interface(pdb_id, chain_id)
    print(f"select foldx_interface, chain B and resi {' + '.join(residue_list)}")


if __name__ == "__main__":
    main()
