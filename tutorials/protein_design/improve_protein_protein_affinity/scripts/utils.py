from typing import Dict, List, Tuple

from Bio.PDB.PDBParser import PDBParser


class AntibodyConstants:
    """Constants related to antibody structure and amino acids."""

    # Three-letter to one-letter amino acid code mapping
    D3TO1: Dict[str, str] = {
        "ALA": "A",
        "CYS": "C",
        "ASP": "D",
        "GLU": "E",
        "PHE": "F",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LYS": "K",
        "LEU": "L",
        "MET": "M",
        "ASN": "N",
        "PRO": "P",
        "GLN": "Q",
        "ARG": "R",
        "SER": "S",
        "THR": "T",
        "VAL": "V",
        "TRP": "W",
        "TYR": "Y",
    }

    # Three-letter amino acid names in alphabetical order
    AA_LIST: List[str] = [
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

    # Complementarity-determining regions (CDRs) position ranges (Kabat numbering)
    CDR1: range = range(26, 34)
    CDR2: range = range(51, 59)
    CDR3: range = range(97, 113)


# Legacy constants for backwards compatibility
D3TO1 = AntibodyConstants.D3TO1
AA_LIST = AntibodyConstants.AA_LIST
CDR1 = AntibodyConstants.CDR1
CDR2 = AntibodyConstants.CDR2
CDR3 = AntibodyConstants.CDR3


def extract_sequence_from_pdb(pdb_file: str, chain_id: str) -> List[Tuple[str, str, str]]:
    """
    Extract amino acid sequence from a PDB file for a specific chain.

    Args:
        pdb_file: Path to the PDB file
        chain_id: Chain identifier (e.g., 'A', 'B')

    Returns:
        List of tuples containing (amino_acid_code, chain_id, position)

    Raises:
        FileNotFoundError: If PDB file doesn't exist
        ValueError: If chain_id is not found in the structure
    """
    if not pdb_file:
        raise ValueError("PDB file path cannot be empty")
    if not chain_id:
        raise ValueError("Chain ID cannot be empty")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(chain_id, pdb_file)

    sequence = []
    chain_found = False

    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                chain_found = True
                for pos, residue in enumerate(chain, 1):
                    if residue.resname in AntibodyConstants.D3TO1:
                        amino_acid = AntibodyConstants.D3TO1[residue.resname]
                        sequence.append((amino_acid, chain.id, str(pos)))

    if not chain_found:
        raise ValueError(f"Chain '{chain_id}' not found in structure")

    return sequence
