#!/usr/bin/env python3
"""
Script to repair PDB structures using FoldX.

Usage:
    python run_foldx_repair.py <pdb_id>

Args:
    pdb_id: PDB identifier (without .pdb extension)
"""

import sys

from foldx_wrapper import FoldXRunner


def main() -> None:
    """Main function to run FoldX RepairPDB command."""
    if len(sys.argv) != 2:
        print("Usage: python run_foldx_repair.py <pdb_id>")
        sys.exit(1)

    pdb_id = sys.argv[1]

    if not pdb_id:
        print("Error: PDB ID cannot be empty")
        sys.exit(1)

    try:
        runner = FoldXRunner()
        result = runner.repair_pdb(pdb_id)
        print(f"Successfully repaired PDB: {pdb_id}")
        if result.stdout:
            print("FoldX output:", result.stdout)
    except Exception as e:
        print(f"Error repairing PDB {pdb_id}: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
