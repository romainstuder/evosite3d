import os
import sys

pdb_id = sys.argv[1]


def main():
    os.system(f"wget --no-clobber https://files.rcsb.org/download/{pdb_id}.pdb")


if __name__ == "__main__":
    main()
