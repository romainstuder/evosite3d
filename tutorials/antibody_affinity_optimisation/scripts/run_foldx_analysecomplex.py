import os
import sys

pdb_id = sys.argv[1]
chain1 = sys.argv[2]
chain2 = sys.argv[3]

# https://foldxsuite.crg.eu/command/AnalyseComplex


def main():
    os.system(
        f"./foldx --command=AnalyseComplex \
                        --pdb={pdb_id}_Repair.pdb \
                        --analyseComplexChains={chain1},{chain2} \
                        --complexWithDNA=false"
    )


if __name__ == "__main__":
    main()
