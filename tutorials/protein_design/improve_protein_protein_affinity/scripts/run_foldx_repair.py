#!/usr/bin/env python3
"""
Script to repair PDB structures using FoldX.

Usage:
    python run_foldx_repair.py <pdb_id>

Args:
    pdb_id: PDB identifier (without .pdb extension)
"""

import argparse
import logging
import sys

from foldx_wrapper import FoldXRunner

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main() -> None:
    """Main function to run FoldX RepairPDB command."""
    parser = argparse.ArgumentParser(description="Repair PDB structures using FoldX")
    parser.add_argument("pdb_id", help="PDB identifier (without .pdb extension)")

    args = parser.parse_args()

    pdb_id = args.pdb_id

    try:
        runner = FoldXRunner()
        result = runner.repair_pdb(pdb_id)
        logger.info(f"Successfully repaired PDB: {pdb_id}")
        if result.stdout:
            logger.debug("FoldX output: %s", result.stdout)
    except Exception as e:
        logger.error(f"Error repairing PDB {pdb_id}: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
