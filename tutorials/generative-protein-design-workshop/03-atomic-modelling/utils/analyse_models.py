# analyze_models.py
import glob

import numpy as np
import pandas as pd
from Bio import PDB


def analyze_model(pdb_file):
    """Analyze a single model"""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_file)

    # Count atoms and residues
    n_atoms = len(list(structure.get_atoms()))
    n_residues = len(list(structure.get_residues()))

    # Check for backbone completeness
    backbone_complete = True
    for residue in structure.get_residues():
        for atom in ["N", "CA", "C", "O"]:
            if atom not in residue:
                backbone_complete = False
                break

    # Calculate radius of gyration
    coords = []
    for atom in structure.get_atoms():
        coords.append(atom.get_coord())
    coords = np.array(coords)
    center = np.mean(coords, axis=0)
    rg = np.sqrt(np.mean(np.sum((coords - center) ** 2, axis=1)))

    return {
        "file": pdb_file,
        "n_atoms": n_atoms,
        "n_residues": n_residues,
        "backbone_complete": backbone_complete,
        "radius_of_gyration": rg,
    }


# Analyze all generated models
print("Analyzing MODELLER output...")
models = glob.glob("target.B*.pdb")
results = []

for model in models:
    results.append(analyze_model(model))

df = pd.DataFrame(results)
print("\nModel Analysis:")
print(df.to_string(index=False))
print(f"\nTotal models generated: {len(models)}")
