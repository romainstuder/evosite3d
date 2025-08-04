#!/usr/bin/env python3

import argparse
import os

from Bio import PDB, SeqIO


def extract_ca_trace(pdb_file, output_file):
    """Extract only CA atoms from structure"""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("scaffold", pdb_file)

    class CASelect(PDB.Select):
        def accept_atom(self, atom):
            return atom.get_name() == "CA"

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_file, CASelect())
    print(f"CA trace saved to {output_file}")


def prepare_alignment_file(sequence_file, template_pdb, output_ali):
    """Create alignment file for MODELLER"""
    # Read the designed sequence
    seq_record = list(SeqIO.parse(sequence_file, "fasta"))[0]
    sequence = str(seq_record.seq)

    # Get template sequence from PDB
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("template", template_pdb)

    template_seq = ""
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[0] == " ":  # Standard residue
                    template_seq += "G"  # Use G for unknown in template

    # Create alignment file
    with open(output_ali, "w") as f:
        f.write(">P1;template_ca\n")
        f.write("structure:template_ca:FIRST:A:LAST:A::::\n")
        f.write(f"{template_seq}*\n\n")

        f.write(">P1;target\n")
        f.write("sequence:target::::::::\n")
        f.write(f"{sequence}*\n")

    print(f"Alignment file created: {output_ali}")


def main():
    """Main function with command line argument parsing"""
    parser = argparse.ArgumentParser(description="Prepare inputs for protein structure modeling")
    parser.add_argument("--scaffold-pdb", required=True, help="Input scaffold PDB file")
    parser.add_argument(
        "--sequence-file", required=True, help="Input FASTA file with designed sequence"
    )
    parser.add_argument(
        "--ca-output",
        default="template_ca.pdb",
        help="Output file for CA trace (default: template_ca.pdb)",
    )
    parser.add_argument(
        "--alignment-output",
        default="alignment.ali",
        help="Output alignment file for MODELLER (default: alignment.ali)",
    )
    parser.add_argument("--verbose", action="store_true", help="Print detailed information")

    args = parser.parse_args()

    # Check if input files exist
    if not os.path.exists(args.scaffold_pdb):
        print(f"Error: Scaffold PDB file {args.scaffold_pdb} not found")
        return

    if not os.path.exists(args.sequence_file):
        print(f"Error: Sequence file {args.sequence_file} not found")
        return

    if args.verbose:
        print(f"Processing scaffold: {args.scaffold_pdb}")
        print(f"Processing sequence: {args.sequence_file}")

    # Prepare files
    extract_ca_trace(args.scaffold_pdb, args.ca_output)
    prepare_alignment_file(args.sequence_file, args.ca_output, args.alignment_output)

    if args.verbose:
        print("Input preparation completed successfully")


if __name__ == "__main__":
    main()
