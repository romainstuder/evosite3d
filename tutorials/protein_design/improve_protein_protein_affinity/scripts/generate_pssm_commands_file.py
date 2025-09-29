#!/usr/bin/env python3
"""
Generate PSSM command files for FoldX batch processing.

This script generates a list of FoldX PSSM commands for each position in a protein chain,
allowing for parallel execution of mutagenesis scans.

Usage:
    python generate_pssm_commands_file.py <pdb_id> <other_chain> <target_chain>

Args:
    pdb_id: PDB identifier (without .pdb extension)
    other_chain: Chain identifier for the binding partner
    target_chain: Chain identifier for the target to mutate

Output:
    command_list.txt: File containing all PSSM commands
"""

import os
import sys
from pathlib import Path
from typing import List, Tuple

from utils import extract_sequence_from_pdb


def create_output_directories(pdb_id: str, sequence: List[Tuple[str, str, str]]) -> None:
    """
    Create output directories for PSSM results.

    Args:
        pdb_id: PDB identifier
        sequence: List of residue information tuples
    """
    base_dir = Path("./pssm_results/")
    base_dir.mkdir(exist_ok=True)

    for residue in sequence:
        position = residue[2].zfill(3)
        output_dir = base_dir / f"{pdb_id}_Repair_{position}"
        output_dir.mkdir(exist_ok=True)


def generate_pssm_commands(
    pdb_id: str, other_chain: str, target_chain: str, sequence: List[Tuple[str, str, str]]
) -> List[str]:
    """
    Generate PSSM command strings for each residue position.

    Args:
        pdb_id: PDB identifier
        other_chain: Chain identifier for binding partner
        target_chain: Chain identifier for target chain
        sequence: List of residue information tuples

    Returns:
        List of command strings
    """
    commands = []

    for residue in sequence:
        res_to_mutate = "".join(residue)
        position = residue[2].zfill(3)
        output_dir = f"./pssm_results/{pdb_id}_Repair_{position}"

        command_parts = [
            "./foldx",
            "--command=Pssm",
            f"--analyseComplexChains={other_chain},{target_chain}",
            f"--pdb={pdb_id}_Repair.pdb",
            f"--positions={res_to_mutate}a",
            f"--output-dir={output_dir}",
        ]

        commands.append(" ".join(command_parts))

    return commands


def write_command_file(commands: List[str], filename: str = "command_list.txt") -> None:
    """
    Write commands to a file.

    Args:
        commands: List of command strings
        filename: Output filename
    """
    with open(filename, "w") as f:
        for command in commands:
            f.write(command + "\n")


def main() -> None:
    """Main function to generate PSSM command file."""
    if len(sys.argv) != 4:
        print("Usage: python generate_pssm_commands_file.py <pdb_id> <other_chain> <target_chain>")
        sys.exit(1)

    pdb_id, other_chain, target_chain = sys.argv[1:4]

    if not all([pdb_id, other_chain, target_chain]):
        print("Error: All arguments must be provided")
        sys.exit(1)

    pdb_file = f"{pdb_id}_Repair.pdb"

    if not os.path.exists(pdb_file):
        print(f"Error: PDB file not found: {pdb_file}")
        sys.exit(1)

    try:
        # Extract sequence information
        sequence = extract_sequence_from_pdb(pdb_file, target_chain)

        if not sequence:
            print(f"Error: No residues found for chain {target_chain}")
            sys.exit(1)

        # Create output directories
        create_output_directories(pdb_id, sequence)

        # Generate commands
        commands = generate_pssm_commands(pdb_id, other_chain, target_chain, sequence)

        # Write command file
        write_command_file(commands)

        print(f"Generated {len(commands)} PSSM commands")
        print("To run in parallel (7 processes):")
        print("cat command_list.txt | xargs -P 7 -I {} sh -c '{}'")

    except Exception as e:
        print(f"Error generating PSSM commands: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
