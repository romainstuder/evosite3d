#!/usr/bin/env python3

from Bio.PDB import PDBIO, PDBParser, Select


class ResidueSelect(Select):
    def __init__(self, chain, start, end):
        self.chain = chain
        self.start = start
        self.end = end

    def accept_residue(self, residue):
        return residue.parent.id == self.chain and self.start <= residue.id[1] <= self.end


def extract_motif(input_pdb, chain, start, end, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, ResidueSelect(chain, start, end))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--chain", default="A")
    parser.add_argument("--residues", required=True, help="start-end format")
    parser.add_argument("--output", required=True)

    args = parser.parse_args()
    start, end = map(int, args.residues.split("-"))
    extract_motif(args.input, args.chain, start, end, args.output)
