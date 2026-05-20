#!/usr/bin/env python3

import numpy as np
from Bio.PDB import PDBIO, PDBParser


def translate_coordinates(input_pdb, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)

    # Find minimum coordinates
    min_coords = np.array([float("inf"), float("inf"), float("inf")])

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coords = atom.get_coord()
                    min_coords = np.minimum(min_coords, coords)

    # Calculate translation vector (add small buffer to ensure positive values)
    translation = -min_coords + 1.0

    # Apply translation
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coords = atom.get_coord()
                    new_coords = coords + translation
                    atom.set_coord(new_coords)

    # Save translated structure
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

    print(f"Translation applied: {translation}")
    print(f"Saved translated structure to: {output_pdb}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input PDB file")
    parser.add_argument("--output", required=True, help="Output PDB file")

    args = parser.parse_args()
    translate_coordinates(args.input, args.output)
