import argparse
import os


def main():
    parser = argparse.ArgumentParser(description="Download PDB file from RCSB")
    parser.add_argument("pdb_id", help="PDB identifier")
    parser.add_argument("workdir", help="Working directory")

    args = parser.parse_args()

    pdb_url = "https://files.rcsb.org/download"
    pdb_file = f"{args.pdb_id}.pdb"
    workdir = args.workdir
    pdb_path = f" {workdir}/{pdb_file}"

    # checking if the workdir directory exists or not.
    if not os.path.exists(workdir):
        # if the workdir directory is not present, then create it.
        os.makedirs(workdir)
    os.system(f"wget --no-clobber --output-document {pdb_path}  {pdb_url}/{pdb_file}")


if __name__ == "__main__":
    main()
