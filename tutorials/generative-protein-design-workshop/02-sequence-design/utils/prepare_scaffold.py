#!/usr/bin/env python3

import argparse

from Bio import PDB


def clean_pdb_for_mpnn(input_pdb, output_pdb):
    """Clean and prepare PDB file for ProteinMPNN"""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("scaffold", input_pdb)

    # Remove heteroatoms and water
    class CleanSelect(PDB.Select):
        def accept_residue(self, residue):
            return residue.get_id()[0] == " "

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, CleanSelect())

    print(f"Cleaned structure saved to {output_pdb}")


def main():
    parser = argparse.ArgumentParser(description="Clean and prepare PDB file for ProteinMPNN")
    parser.add_argument("--input_pdb", help="Input PDB file path")
    parser.add_argument("--output_pdb", help="Output cleaned PDB file path")

    args = parser.parse_args()
    clean_pdb_for_mpnn(args.input_pdb, args.output_pdb)


if __name__ == "__main__":
    main()
