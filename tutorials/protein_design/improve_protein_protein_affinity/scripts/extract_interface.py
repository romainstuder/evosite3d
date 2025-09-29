import argparse
import logging

import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


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
    parser = argparse.ArgumentParser(description="Extract interface residues from FoldX output")
    parser.add_argument("pdb_id", help="PDB identifier (without .pdb extension)")
    parser.add_argument("chain_id", help="Chain identifier")

    args = parser.parse_args()

    pdb_id = args.pdb_id
    chain_id = args.chain_id

    residue_list = extract_interface(pdb_id, chain_id)
    logger.info(f"select foldx_interface, chain B and resi {' + '.join(residue_list)}")


if __name__ == "__main__":
    main()
