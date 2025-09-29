import argparse
import os


def main():
    parser = argparse.ArgumentParser(description="Download PDB file from RCSB")
    parser.add_argument("pdb_id", help="PDB identifier")

    args = parser.parse_args()

    pdb_id = args.pdb_id

    os.system(f"wget --no-clobber https://files.rcsb.org/download/{pdb_id}.pdb")


if __name__ == "__main__":
    main()
