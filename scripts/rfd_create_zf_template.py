#!/usr/bin/env python3

from Bio.PDB import PDBIO, Atom, Chain, Model, Residue, Structure


def create_zf_template():
    """Create a minimal zinc-finger template structure"""

    # Realistic C2H2 zinc finger coordinates based on typical zinc finger geometry
    # Zinc ion positioned at origin with tetrahedral coordination
    zf_coords = {
        "CYS_1": {"CA": [-2.1, 1.8, -1.2], "CB": [-1.2, 0.9, -2.1], "SG": [0.3, 1.4, -1.8]},
        "CYS_2": {"CA": [1.8, -2.1, -1.2], "CB": [0.9, -1.2, -2.1], "SG": [-0.3, -1.4, -1.8]},
        "HIS_1": {
            "CA": [-1.8, -1.8, 1.8],
            "CB": [-0.9, -0.9, 2.7],
            "CG": [0.2, -1.2, 2.1],
            "ND1": [0.8, -0.5, 1.2],
            "NE2": [1.1, -2.1, 2.4],
        },
        "HIS_2": {
            "CA": [2.1, 1.8, 1.2],
            "CB": [1.2, 0.9, 2.1],
            "CG": [-0.2, 1.2, 2.1],
            "ND1": [-0.8, 0.5, 1.2],
            "NE2": [-1.1, 2.1, 2.4],
        },
    }

    zf_structure = Structure.Structure("ZF_template")
    model = Model.Model(0)
    chain = Chain.Chain("A")

    for i, (res_name, coords) in enumerate(zf_coords.items(), 1):
        residue = Residue.Residue((" ", i, " "), res_name[:3], " ")
        for atom_name, coord in coords.items():
            atom = Atom.Atom(atom_name, coord, 0.0, 1.0, " ", atom_name, i)
            residue.add(atom)
        chain.add(residue)

    model.add(chain)
    zf_structure.add(model)
    return zf_structure


# Save template
structure = create_zf_template()
io = PDBIO()
io.set_structure(structure)
io.save("zf_template.pdb")
