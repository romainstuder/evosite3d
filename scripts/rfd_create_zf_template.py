#!/usr/bin/env python3

from Bio.PDB import PDBIO, Atom, Chain, Model, Residue, Structure


def create_zf_template():
    """Create a minimal zinc-finger template structure"""

    # Idealized C2H2 zinc finger coordinates (simplified)
    zf_coords = {
        "CYS_1": {"CA": [0.0, 0.0, 0.0], "CB": [1.5, 0.0, 0.0], "SG": [2.0, 1.5, 0.0]},
        "CYS_2": {"CA": [8.0, 2.0, 1.0], "CB": [9.0, 2.5, 1.0], "SG": [10.0, 1.8, 2.0]},
        "HIS_1": {"CA": [15.0, 8.0, 5.0], "CB": [15.5, 9.0, 6.0], "CG": [16.0, 8.5, 7.0]},
        "HIS_2": {"CA": [18.0, 12.0, 8.0], "CB": [18.5, 13.0, 7.0], "CG": [19.0, 12.5, 6.0]},
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
