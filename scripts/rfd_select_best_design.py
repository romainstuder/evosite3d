#!/usr/bin/env python3

from Bio import PDB
from rfd_analyse_structures import calculate_helix_content, calculate_radius_of_gyration


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


def evaluate_design(pdb_file):
    """Score a design based on multiple criteria"""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    # Calculate metrics
    scores = {
        "helix_content": calculate_helix_content(pdb_file),
        "compactness": calculate_radius_of_gyration(structure),
        "clash_score": check_clashes(structure),
    }

    # Weighted score
    total_score = (
        scores["helix_content"] * 0.4
        + (50 - scores["compactness"]) * 0.3
        + (100 - scores["clash_score"]) * 0.3
    )

    return total_score, scores


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

# Evaluate all designs and select best
best_design = max(designs, key=lambda x: evaluate_design(x["file"])[0])
print(f"Best design: {best_design['file']}")
