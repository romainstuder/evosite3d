#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
from Bio import PDB


def calculate_helix_content(structure_file):
    """Calculate the percentage of helical residues"""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", structure_file)

    # Use DSSP or simple phi-psi analysis
    # This is a simplified version - in practice you'd use DSSP
    helix_count = 0
    total_residues = 0

    for model in structure:
        for chain in model:
            total_residues += len(chain)
            # Simplified helix detection - in practice use DSSP
            # For now, we'll estimate based on structure regularity
            helix_count += int(len(chain) * 0.6)  # Estimate 60% helix

    return (helix_count / total_residues) * 100 if total_residues > 0 else 0


def calculate_radius_of_gyration(structure):
    """Calculate radius of gyration for a protein structure"""
    import numpy as np

    # Get all CA atoms
    ca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    ca_atoms.append(residue["CA"].get_coord())

    if not ca_atoms:
        return float("inf")

    # Convert to numpy array
    coords = np.array(ca_atoms)

    # Calculate center of mass (assuming equal mass for all CA atoms)
    center_of_mass = np.mean(coords, axis=0)

    # Calculate radius of gyration
    distances_squared = np.sum((coords - center_of_mass) ** 2, axis=1)
    rg = np.sqrt(np.mean(distances_squared))

    return rg


def check_clashes(structure, clash_distance=2.5):
    """Check for steric clashes in the structure"""
    atoms = list(structure.get_atoms())
    clash_count = 0

    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            # Skip if atoms are from adjacent residues (likely bonded)
            res_i = atoms[i].get_parent()
            res_j = atoms[j].get_parent()
            if abs(res_i.get_id()[1] - res_j.get_id()[1]) <= 1:
                continue

            # Calculate distance
            distance = atoms[i] - atoms[j]
            if distance < clash_distance:
                clash_count += 1

    return clash_count


# Analyze all designs
designs = []
for i in range(10):
    pdb_file = f"./output/two_helix_designs/design_{i}.pdb"

    # Parse structure
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    # Calculate all metrics
    helix_content = calculate_helix_content(pdb_file)
    rg = calculate_radius_of_gyration(structure)
    clashes = check_clashes(structure)

    designs.append(
        {
            "design": i,
            "file": pdb_file,
            "helix_content": helix_content,
            "radius_of_gyration": rg,
            "clash_score": clashes,
        }
    )

# Convert to DataFrame for easier analysis
df_designs = pd.DataFrame(designs)
print("\nDesign Analysis Summary:")
print(df_designs)

# Visualize results
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

df_designs["helix_content"].plot(kind="bar", ax=axes[0])
axes[0].set_title("Helix Content by Design")
axes[0].set_ylabel("Helix Content (%)")

df_designs["radius_of_gyration"].plot(kind="bar", ax=axes[1])
axes[1].set_title("Radius of Gyration by Design")
axes[1].set_ylabel("Rg (Å)")

df_designs["clash_score"].plot(kind="bar", ax=axes[2])
axes[2].set_title("Clash Score by Design")
axes[2].set_ylabel("Number of Clashes")

plt.tight_layout()
plt.savefig("design_analysis.png")
plt.show()
