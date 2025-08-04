#!/usr/bin/env python3

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


# Clean your scaffold
clean_pdb_for_mpnn("two_helix_scaffold.pdb", "scaffold_clean.pdb")
