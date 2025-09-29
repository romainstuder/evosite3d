#!/usr/bin/env python3
"""
Script to analyze protein complexes using FoldX.

Usage:
    python run_foldx_analysecomplex.py <pdb_id> <chain1> <chain2>

Args:
    pdb_id: PDB identifier (without .pdb extension)
    chain1: First chain identifier
    chain2: Second chain identifier

Reference:
    https://foldxsuite.crg.eu/command/AnalyseComplex
"""

import sys

from foldx_wrapper import FoldXRunner


def main() -> None:
    """Main function to run FoldX AnalyseComplex command."""
    if len(sys.argv) != 4:
        print("Usage: python run_foldx_analysecomplex.py <pdb_id> <chain1> <chain2>")
        sys.exit(1)

    pdb_id, chain1, chain2 = sys.argv[1:4]

    if not all([pdb_id, chain1, chain2]):
        print("Error: All arguments (pdb_id, chain1, chain2) must be provided")
        sys.exit(1)

    try:
        runner = FoldXRunner()
        result = runner.analyse_complex(pdb_id, chain1, chain2)
        print(f"Successfully analyzed complex: {pdb_id} (chains {chain1}, {chain2})")
        if result.stdout:
            print("FoldX output:", result.stdout)
    except Exception as e:
        print(f"Error analyzing complex {pdb_id}: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
