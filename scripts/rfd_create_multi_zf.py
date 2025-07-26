#!/usr/bin/env python3

import argparse

import numpy as np
from Bio.PDB import PDBIO, Atom, Chain, Model, Residue, Structure


def create_multi_zf(num_fingers, finger_spacing, output_file):
    structure = Structure.Structure("Multi_ZF")
    model = Model.Model(0)
    chain = Chain.Chain("A")

    # Standard C2H2 zinc finger template (30 residues)
    zf_template = "PYKCELCGKSFRSDSSLTKHQRIHTGEKPFA"

    residue_id = 1
    x_offset = 0

    for finger in range(num_fingers):
        for i, aa in enumerate(zf_template.strip()):
            # Simple helical geometry for demonstration
            x = x_offset + i * 1.5
            y = 5 * np.sin(i * 0.3)
            z = 5 * np.cos(i * 0.3)

            residue = Residue.Residue((" ", residue_id, " "), aa, " ")
            atom = Atom.Atom("CA", [x, y, z], 0.0, 1.0, " ", "CA", residue_id)
            residue.add(atom)
            chain.add(residue)
            residue_id += 1

        x_offset += 45  # Space between fingers

        # Add linker residues
        if finger < num_fingers - 1:
            linker = "TGEKPFA"  # 7-residue linker
            for i, aa in enumerate(linker):
                x = x_offset + i * 1.5
                y = 0
                z = 0

                residue = Residue.Residue((" ", residue_id, " "), aa, " ")
                atom = Atom.Atom("CA", [x, y, z], 0.0, 1.0, " ", "CA", residue_id)
                residue.add(atom)
                chain.add(residue)
                residue_id += 1

    model.add(chain)
    structure.add(model)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fingers", type=int, default=3)
    parser.add_argument("--spacing", type=int, default=7)
    parser.add_argument("--output", default="multi_zf.pdb")
    args = parser.parse_args()

    create_multi_zf(args.fingers, args.spacing, args.output)
