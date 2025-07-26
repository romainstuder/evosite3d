#!/usr/bin/env python3

import numpy as np
from Bio.PDB import PDBParser


def analyze_zf_design(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("design", pdb_file)

    # Find zinc-coordinating residues
    zn_coords = []
    for residue in structure.get_residues():
        if residue.resname in ["CYS", "HIS"]:
            if residue.resname == "CYS" and "SG" in residue:
                zn_coords.append(residue["SG"].coord)
            elif residue.resname == "HIS" and "NE2" in residue:
                zn_coords.append(residue["NE2"].coord)

    if len(zn_coords) >= 4:
        # Calculate potential zinc position (centroid)
        zn_pos = np.mean(zn_coords[:4], axis=0)

        # Calculate distances to coordinating atoms
        distances = [np.linalg.norm(coord - zn_pos) for coord in zn_coords[:4]]

        print(f"Design: {pdb_file}")
        print(f"Zinc coordination distances: {distances}")
        print(f"Average distance: {np.mean(distances):.2f} Å")
        print(f"Distance std: {np.std(distances):.2f} Å")
        print("Good zinc coordination: 2.0-2.5 Å average")
        print()


# Analyze all designs
for i in range(5):
    analyze_zf_design(f"zf_scaffold_example_{i}.pdb")
