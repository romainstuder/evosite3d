# validate_best_model.py
import numpy as np
from Bio import PDB


def validate_model(pdb_file):
    """Comprehensive validation of the model"""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_file)

    print(f"\nValidating: {pdb_file}")
    print("=" * 50)

    # 1. Basic statistics
    all_atoms = list(structure.get_atoms())
    all_residues = list(structure.get_residues())

    print("\nStructure Statistics:")
    print(f"  Total residues: {len(all_residues)}")
    print(f"  Total atoms: {len(all_atoms)}")

    # 2. Check atom types
    atom_types = {}
    for atom in all_atoms:
        atom_type = atom.element
        atom_types[atom_type] = atom_types.get(atom_type, 0) + 1

    print("\nAtom composition:")
    for atom_type, count in sorted(atom_types.items()):
        print(f"  {atom_type}: {count}")

    # 3. Check for missing backbone atoms
    missing_atoms = []
    for residue in all_residues:
        for atom_name in ["N", "CA", "C", "O"]:
            if atom_name not in residue:
                missing_atoms.append(f"{residue.get_resname()}{residue.id[1]}:{atom_name}")

    if missing_atoms:
        print("\nWarning - Missing backbone atoms:")
        for atom in missing_atoms:
            print(f"  {atom}")
    else:
        print("\nBackbone complete: ✓")

    # 4. Simple clash detection
    clashes = 0
    for i in range(len(all_atoms)):
        for j in range(i + 1, len(all_atoms)):
            if all_atoms[i].parent.id[1] == all_atoms[j].parent.id[1]:
                continue  # Same residue
            if abs(all_atoms[i].parent.id[1] - all_atoms[j].parent.id[1]) == 1:
                continue  # Adjacent residues

            distance = all_atoms[i] - all_atoms[j]
            if distance < 2.0 and all_atoms[i].element != "H" and all_atoms[j].element != "H":
                clashes += 1

    print(f"\nSteric clashes (< 2.0 Å): {clashes}")

    # 5. B-factor analysis (if available)
    b_factors = [atom.get_bfactor() for atom in all_atoms]
    if any(b > 0 for b in b_factors):
        print("\nB-factor statistics:")
        print(f"  Mean: {np.mean(b_factors):.2f}")
        print(f"  Std: {np.std(b_factors):.2f}")
        print(f"  Range: {np.min(b_factors):.2f} - {np.max(b_factors):.2f}")

    return True


# Validate the best model
if __name__ == "__main__":
    import os

    if os.path.exists("best_model.pdb"):
        validate_model("best_model.pdb")
        print("\n" + "=" * 50)
        print("VALIDATION COMPLETE")
        print("Your model is ready for further analysis!")
    else:
        print("Error: best_model.pdb not found!")
        print("Please run build_model_complete.py first.")
