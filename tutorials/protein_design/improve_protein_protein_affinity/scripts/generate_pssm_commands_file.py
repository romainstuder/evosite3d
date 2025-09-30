#!/usr/bin/env python3
"""Generate PSSM command files for FoldX batch processing.

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

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Tuple

import pandas as pd
from utils import extract_sequence_from_pdb

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def extract_interface(pdb_id: str, chain_id: str, workdir: str) -> List[str]:
    """
    Extract interface residue positions from FoldX AnalyseComplex output.

    Args:
        pdb_id: PDB identifier
        chain_id: Chain identifier
        workdir: Working directory

    Returns:
        List of residue positions as strings (e.g., ['45', '47', '52'])
    """
    workdir_path = Path(workdir)
    interface_file = workdir_path / f"Interface_Residues_{pdb_id}_Repair_AC.fxout"

    if not interface_file.exists():
        logger.error(f"Interface file not found: {interface_file}")
        logger.error("Please run run_foldx_analysecomplex.py first")
        sys.exit(1)

    df = pd.read_csv(interface_file, sep="\t", skiprows=10, header=0)
    interface_residues = df.columns.tolist()
    interface_residues.remove(df.columns[-1])

    pos_in_chain = {}

    for res in interface_residues:
        chain = res[1]
        pos = res[2:]
        pos_in_chain.setdefault(chain, []).append(pos)

    if chain_id not in pos_in_chain:
        logger.error(f"No interface residues found for chain {chain_id}")
        logger.error(f"Available chains: {list(pos_in_chain.keys())}")
        sys.exit(1)

    return pos_in_chain[chain_id]


def create_output_directories(
    pdb_id: str, sequence: List[Tuple[str, str, str]], workdir: str
) -> None:
    """
    Create output directories for PSSM results.

    Args:
        pdb_id: PDB identifier
        sequence: List of residue information tuples
        workdir: Working directory
    """
    base_dir = Path(workdir) / "pssm_results"
    base_dir.mkdir(exist_ok=True)

    for residue in sequence:
        position = residue[2].zfill(3)
        output_dir = base_dir / f"{pdb_id}_Repair_{position}"
        output_dir.mkdir(exist_ok=True)


def generate_pssm_commands(
    pdb_id: str,
    other_chain: str,
    target_chain: str,
    sequence: List[Tuple[str, str, str]],
    workdir: str,
) -> List[str]:
    """
    Generate PSSM command strings for each residue position.

    Args:
        pdb_id: PDB identifier
        other_chain: Chain identifier for binding partner
        target_chain: Chain identifier for target chain
        sequence: List of residue information tuples
        workdir: Working directory

    Returns:
        List of command strings
    """
    commands = []
    workdir_path = Path(workdir)

    for residue in sequence:
        res_to_mutate = "".join(residue)
        position = residue[2].zfill(3)
        output_dir = workdir_path / "pssm_results" / f"{pdb_id}_Repair_{position}"

        command_parts = [
            "./foldx",
            "--command=Pssm",
            f"--analyseComplexChains={other_chain},{target_chain}",
            f"--pdb={pdb_id}_Repair.pdb",
            f"--pdb-dir={workdir}",
            f"--positions={res_to_mutate}a",
            f"--output-dir={output_dir}",
        ]

        commands.append(" ".join(command_parts))

    return commands


def write_command_file(
    commands: List[str], workdir: str, filename: str = "command_list.txt"
) -> None:
    """
    Write commands to a file.

    Args:
        commands: List of command strings
        workdir: Working directory
        filename: Output filename
    """
    output_file = Path(workdir) / filename
    output_file.write_text("\n".join(commands) + "\n")


def main() -> None:
    """Main function to generate PSSM command file."""
    parser = argparse.ArgumentParser(
        description="Generate PSSM command files for FoldX batch processing"
    )
    parser.add_argument("pdb_id", help="PDB identifier (without .pdb extension)")
    parser.add_argument("other_chain", help="Chain identifier for the binding partner")
    parser.add_argument("target_chain", help="Chain identifier for the target to mutate")
    parser.add_argument(
        "--full-sequence",
        action="store_true",
        help="Scan all residues in the chain (default: interface residues only)",
    )
    parser.add_argument("workdir", help="Working directory")

    args = parser.parse_args()

    pdb_id = args.pdb_id
    other_chain = args.other_chain
    target_chain = args.target_chain
    workdir = args.workdir

    workdir_path = Path(workdir)
    pdb_file = workdir_path / f"{pdb_id}_Repair.pdb"

    if not pdb_file.exists():
        logger.error(f"Error: PDB file not found: {pdb_file}")
        sys.exit(1)

    try:
        # Extract sequence information
        sequence = extract_sequence_from_pdb(str(pdb_file), target_chain)

        if not sequence:
            logger.error(f"Error: No residues found for chain {target_chain}")
            sys.exit(1)

        # Filter to interface residues by default (unless --full-sequence is specified)
        if not args.full_sequence:
            logger.info(
                "Using interface residues only (use --full-sequence to scan all residues)..."
            )
            interface_positions = extract_interface(pdb_id, target_chain, workdir)
            interface_positions_set = set(interface_positions)

            # Filter sequence to only include interface positions
            filtered_sequence = [
                residue for residue in sequence if residue[2] in interface_positions_set
            ]

            if not filtered_sequence:
                logger.error("No interface residues found in sequence")
                sys.exit(1)

            logger.info(
                f"Filtered from {len(sequence)} to {len(filtered_sequence)} interface residues"
            )
            logger.info(f"Interface positions: {', '.join(sorted(interface_positions, key=int))}")
            sequence = filtered_sequence
        else:
            logger.info(f"Using all {len(sequence)} residues in chain {target_chain}")

        # Create output directories
        create_output_directories(pdb_id, sequence, workdir)

        # Generate commands
        commands = generate_pssm_commands(pdb_id, other_chain, target_chain, sequence, workdir)

        # Write command file
        command_list_path = workdir_path / "command_list.txt"
        write_command_file(commands, workdir)

        logger.info(f"Generated {len(commands)} PSSM commands")
        logger.info("To run in parallel (7 processes):")
        logger.info(f"cat {command_list_path} | xargs -P 7 -I {{}} sh -c '{{}}'")

    except Exception as e:
        logger.error(f"Error generating PSSM commands: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
