import os
import sys

pdb_id = sys.argv[1]


def main():
    os.system(
        f"./foldx --command=RepairPDB \
            --pdb={pdb_id}.pdb \
            --repair_Interface=ALL"
    )


if __name__ == "__main__":
    main()
