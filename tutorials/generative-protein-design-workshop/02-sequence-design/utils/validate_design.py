#!/usr/bin/env python3

import argparse

import numpy as np
from Bio import PDB, SeqIO
from Bio.PDB import DSSP


def check_sequence_structure_compatibility(
    sequence, dssp, hydrophobic_residues="ILMFWYV", helix_offsets=None
):
    """Check sequence-structure compatibility with configurable parameters"""
    if helix_offsets is None:
        helix_offsets = [0, 3, 4, 7]

    # Check helix positions match hydrophobic pattern
    helix_positions = []
    for i, ss in enumerate(dssp):
        if ss[2] == "H":  # Alpha helix
            helix_positions.append(i)

    # Analyze hydrophobic periodicity in helices
    scores = []
    for i in helix_positions:
        if i < len(sequence):
            # Check helix wheel pattern
            hydrophobic_score = 0
            for offset in helix_offsets:
                if i + offset < len(sequence):
                    if sequence[i + offset] in hydrophobic_residues:
                        hydrophobic_score += 1
            scores.append(hydrophobic_score)

    return np.mean(scores) if scores else 0


def main():
    """Main function to validate designs with command line arguments"""
    parser = argparse.ArgumentParser(
        description="Validate sequence-structure compatibility of protein designs"
    )
    parser.add_argument("--sequences-file", help="Input FASTA file containing sequences")
    parser.add_argument("--structure-file", help="Input PDB file containing structure scaffold")
    parser.add_argument(
        "--hydrophobic-residues",
        default="ILMFWYV",
        help="Hydrophobic residues to check (default: ILMFWYV)",
    )
    parser.add_argument(
        "--helix-offsets",
        nargs="+",
        type=int,
        default=[0, 3, 4, 7],
        help="Helix wheel offsets to check (default: 0 3 4 7)",
    )
    parser.add_argument("--verbose", action="store_true", help="Print detailed analysis")

    args = parser.parse_args()

    # Parse structure
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", args.structure_file)

    # Get secondary structure
    model = structure[0]
    dssp = DSSP(model, args.structure_file)

    # Validate top sequences
    for seq_record in SeqIO.parse(args.sequences_file, "fasta"):
        score = check_sequence_structure_compatibility(
            str(seq_record.seq), dssp, args.hydrophobic_residues, args.helix_offsets
        )
        print(f"{seq_record.id}: Compatibility score = {score:.2f}")

        if args.verbose:
            hydrophobic_content = sum(
                seq_record.seq.count(aa) for aa in args.hydrophobic_residues
            ) / len(seq_record.seq)
            print(f"  Sequence length: {len(seq_record.seq)}")
            print(f"  Hydrophobic content: " f"{hydrophobic_content:.2%}")


if __name__ == "__main__":
    main()
