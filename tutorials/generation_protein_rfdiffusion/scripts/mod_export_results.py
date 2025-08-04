import os
import shutil

from Bio import PDB, SeqIO


def create_project_summary():
    """Create a summary of the entire design project"""

    # Create results directory
    os.makedirs("final_results", exist_ok=True)

    # Copy key files
    files_to_copy = [
        ("two_helix_scaffold.pdb", "Stage 1 - RFdiffusion backbone"),
        ("designed_sequence.fasta", "Stage 2 - ProteinMPNN sequence"),
        ("best_model.pdb", "Stage 3 - MODELLER full model"),
    ]

    for filename, description in files_to_copy:
        if os.path.exists(filename):
            shutil.copy(filename, f"final_results/{filename}")
            print(f"✓ Copied {filename} - {description}")
        else:
            print(f"✗ Missing {filename}")

    # Create summary report
    with open("final_results/design_summary.txt", "w") as f:
        f.write("PROTEIN DESIGN PROJECT SUMMARY\n")
        f.write("=" * 50 + "\n\n")

        # Get sequence
        if os.path.exists("designed_sequence.fasta"):
            seq_record = list(SeqIO.parse("designed_sequence.fasta", "fasta"))[0]
            sequence = str(seq_record.seq)
            f.write(f"Designed Sequence ({len(sequence)} aa):\n")
            f.write(f"{sequence}\n\n")

        # Get structure stats
        if os.path.exists("best_model.pdb"):
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure("model", "best_model.pdb")
            n_atoms = len(list(structure.get_atoms()))
            n_residues = len(list(structure.get_residues()))

            f.write("Final Model Statistics:\n")
            f.write(f"- Residues: {n_residues}\n")
            f.write(f"- Atoms: {n_atoms}\n")

        f.write("\nPipeline:\n")
        f.write("1. RFdiffusion → Backbone generation\n")
        f.write("2. ProteinMPNN → Sequence design\n")
        f.write("3. MODELLER → Full atomic model\n")

        f.write("\nNext Steps:\n")
        f.write("- Validate with AlphaFold2\n")
        f.write("- Run MD simulations\n")
        f.write("- Order gene synthesis\n")
        f.write("- Express and characterize\n")

    print("\n" + "=" * 50)
    print("PROJECT COMPLETE!")
    print("=" * 50)
    print("All results saved in: final_results/")
    print("\nYour designed protein is ready for experimental validation!")


# Run the summary
if __name__ == "__main__":
    create_project_summary()
