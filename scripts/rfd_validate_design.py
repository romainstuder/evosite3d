#!/usr/bin/env python3

import numpy as np
from Bio import PDB, SeqIO
from Bio.PDB import DSSP


def check_sequence_structure_compatibility(sequence, structure_file):
    """Basic checks for sequence-structure compatibility"""

    # Parse structure
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", structure_file)

    # Get secondary structure
    model = structure[0]
    dssp = DSSP(model, structure_file)

    # Check helix positions match hydrophobic pattern
    helix_positions = []
    for i, (res, ss) in enumerate(dssp):
        if ss[2] == "H":  # Alpha helix
            helix_positions.append(i)

    # Analyze hydrophobic periodicity in helices
    scores = []
    for i in helix_positions:
        if i < len(sequence):
            # Check i, i+3, i+4 pattern (helix wheel)
            hydrophobic_score = 0
            for offset in [0, 3, 4, 7]:
                if i + offset < len(sequence):
                    if sequence[i + offset] in "ILMFWYV":
                        hydrophobic_score += 1
            scores.append(hydrophobic_score)

    return np.mean(scores) if scores else 0


# Validate top sequences
for seq_record in SeqIO.parse("filtered_sequences.fasta", "fasta"):
    score = check_sequence_structure_compatibility(str(seq_record.seq), "scaffold_clean.pdb")
    print(f"{seq_record.id}: Compatibility score = {score:.2f}")
