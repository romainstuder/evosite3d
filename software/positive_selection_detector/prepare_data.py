#!/usr/bin/env python3
"""Acquire and prepare all data for CodeML site-model analysis.

Convenience wrapper that runs both sub-steps:
    1a. ``fetch_alignment.py`` — CDS alignment + pruned gene tree
    1b. ``fetch_structure.py`` — 3D structure (AlphaFold or RCSB PDB)

When using Nextflow, these are run as separate parallel processes.
This script is provided for standalone (non-Nextflow) usage.

Example:
    prepare_data.py --gene-symbol HLA-DQB1
    prepare_data.py --gene-symbol HLA-DQB1 --outdir HLA_DQB1 --pdb 1UVQ
"""

from __future__ import annotations

import argparse

from _common import add_common_args  # type: ignore
from ensembl_msa_tree_extract import run_extract_for_symbol as run_fetch_alignment  # type: ignore
from fetch_structure import run_fetch_structure  # type: ignore


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Acquire and prepare data (alignment, tree, "
        "3D structure) for CodeML site-model analysis.",
    )
    add_common_args(ap)
    ap.add_argument(
        "--taxon",
        type=int,
        default=9347,
        help="NCBI taxon ID for the Compara gene-tree prune "
        "(default 9347 = Boreoeutheria). Pass 7742 for Vertebrata.",
    )
    ap.add_argument(
        "--pdb",
        default=None,
        help="Override AlphaFold and fetch this RCSB PDB accession instead.",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    run_fetch_alignment(
        gene_symbol=args.gene_symbol,
        outdir=args.outdir,
        taxon=args.taxon,
        uniprot=args.uniprot,
    )
    run_fetch_structure(
        gene_symbol=args.gene_symbol,
        outdir=args.outdir,
        uniprot=args.uniprot,
        pdb_id=args.pdb,
    )


if __name__ == "__main__":
    main()
