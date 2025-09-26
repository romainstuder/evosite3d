from Bio.PDB.PDBParser import PDBParser

D3TO1 = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
}


AA_LIST = [
    "Ala",
    "Cys",
    "Asp",
    "Glu",
    "Phe",
    "Gly",
    "His",
    "Ile",
    "Lys",
    "Leu",
    "Met",
    "Asn",
    "Pro",
    "Gln",
    "Arg",
    "Ser",
    "Thr",
    "Val",
    "Trp",
    "Tyr",
]


CDR1 = range(26, 34)
CDR2 = range(51, 59)
CDR3 = range(97, 113)


def extract_sequence_from_pdb(pdb_file, chain_id):
    # iterate each model, chain, and residue
    # printing out the sequence for each chain

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(chain_id, pdb_file)

    seq = []
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for pos, residue in enumerate(chain):
                    if residue.resname in D3TO1:
                        seq.append((D3TO1[residue.resname], chain.id, str(pos + 1)))
                # print('>some_header\n',''.join(seq))
    return seq
