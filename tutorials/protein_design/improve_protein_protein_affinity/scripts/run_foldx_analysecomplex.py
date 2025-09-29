#!/usr/bin/env python3
"""
Script to analyze protein complexes using FoldX.

Usage:
    python run_foldx_analysecomplex.py <pdb_id> <chain1> <chain2>

Args:
    pdb_id: PDB identifier (without .pdb extension)
    chain1: First chain identifier
    chain2: Second chain identifier

Reference:
    https://foldxsuite.crg.eu/command/AnalyseComplex
"""

import argparse
import logging
import sys

from foldx_wrapper import FoldXRunner

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main() -> None:
    """Main function to run FoldX AnalyseComplex command."""
    parser = argparse.ArgumentParser(description="Analyze protein complexes using FoldX")
    parser.add_argument("pdb_id", help="PDB identifier (without .pdb extension)")
    parser.add_argument("chain1", help="First chain identifier")
    parser.add_argument("chain2", help="Second chain identifier")

    args = parser.parse_args()

    pdb_id = args.pdb_id
    chain1 = args.chain1
    chain2 = args.chain2

    try:
        runner = FoldXRunner()
        result = runner.analyse_complex(pdb_id, chain1, chain2)
        logger.info(f"Successfully analyzed complex: {pdb_id} (chains {chain1}, {chain2})")
        if result.stdout:
            logger.debug("FoldX output: %s", result.stdout)
    except Exception as e:
        logger.error(f"Error analyzing complex {pdb_id}: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
