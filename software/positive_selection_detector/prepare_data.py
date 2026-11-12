#!/usr/bin/env python3
"""Acquire and prepare all data for CodeML site-model analysis.

Convenience wrapper that runs the three data-preparation sub-steps and
collects every output into a single working directory so the standalone
Step 2 / Step 3 scripts (which read from that one directory) find their
inputs:

    1a. ``ensembl_msa_tree_extract.py``   — Bronze: raw Ensembl tree + MSA + CDS
    1b. ``fetch_structure.py``            — 3D structure (AlphaFold or RCSB PDB)
    2a. ``ensembl_msa_tree_transform.py`` — Silver: cleaned, TrimAl-ed,
                                            PHYLIP-ready alignment + pruned tree

When using Nextflow, these run as separate processes writing into the
Bronze/Silver medallion layout. This wrapper instead writes everything
into one directory (``--outdir``, defaulting to a sanitised form of the
gene symbol) — the layout the standalone CodeML/HyPhy/analyse steps expect.

Example:
    prepare_data.py --gene-symbol HLA-DQB1
    prepare_data.py --gene-symbol HLA-DQB1 --outdir HLA_DQB1 --pdb 1UVQ
"""

from __future__ import annotations

import argparse

from _common import add_common_args, resolve_gene
from ensembl_msa_tree_extract import run_extract_for_symbol as run_fetch_alignment
from ensembl_msa_tree_transform import run_transform_for_symbol
from fetch_structure import run_fetch_structure


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
        "(default 9347 = Eutheria, placental mammals). Pass 7742 for Vertebrata.",
    )
    ap.add_argument(
        "--pdb",
        default=None,
        help="Override AlphaFold and fetch this RCSB PDB accession instead.",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()

    # Resolve once to pin a single working directory (``--outdir`` or, by
    # default, ``./<prefix>``). Every sub-step is pointed at it explicitly so
    # Bronze, structure, and Silver outputs land together — otherwise each
    # step would fall back to a different default root (results/1_bronze,
    # ./<prefix>, results/2_silver) and the standalone Step 2 would not find
    # the alignment.
    ids = resolve_gene(args.gene_symbol, outdir=args.outdir, uniprot_override=args.uniprot)
    workdir = ids.workdir

    # 1a. Bronze — raw Ensembl tree + protein MSA + CDS.
    run_fetch_alignment(
        gene_symbol=args.gene_symbol,
        outdir=workdir,
        taxon=args.taxon,
        uniprot=args.uniprot,
    )
    # 1b. 3D structure. Run before the transform so that an over-gappy
    # alignment (which makes the transform stop early) still leaves the
    # fetched structure on disk.
    run_fetch_structure(
        gene_symbol=args.gene_symbol,
        outdir=workdir,
        uniprot=args.uniprot,
        pdb_id=args.pdb,
    )
    # 2a. Silver — clean/TrimAl/PHYLIP-ready alignment + pruned tree, read
    # from and written back into the same working directory.
    run_transform_for_symbol(
        gene_symbol=args.gene_symbol,
        outdir=workdir,
        indir=workdir,
        uniprot=args.uniprot,
    )


if __name__ == "__main__":
    main()
