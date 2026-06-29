#!/usr/bin/env python3
"""Step 2b — Renumber Bronze PDBs to UniProt residue numbering (Silver).

Reads the ``{prefix}_structures.json`` written by ``fetch_structure.py``
and, for each RCSB PDB listed, applies SIFTS ``pdb_resnum → uniprot_pos``
renumbering on the UniProt-mapped chain. AlphaFold models and PDBs with
no SIFTS data are copied unchanged.

Outputs are written to ``--outdir`` (typically ``2_silver``) using the
same filenames as Bronze, plus a copy of ``{prefix}_structures.json``.

Example:
    renumber_structures.py --gene-symbol HLA-DQB1 \\
        --indir results/1_bronze/HLA_DQB1 \\
        --outdir results/2_silver/HLA_DQB1
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

# ``sifts`` (lib/) is importable because importing ``_common`` puts that
# directory on sys.path.
from _common import add_common_args, log, resolve_gene  # type: ignore
from sifts import (
    build_residue_map,
    renumber_pdb_chain,
)


def renumber_one(
    bronze_pdb: Path,
    silver_pdb: Path,
    uniprot: str,
    pdb_id: str,
    chain_id: str,
) -> None:
    """Renumber a single Bronze PDB into Silver via SIFTS.

    On any error or empty SIFTS map, copies the file unchanged.
    """
    silver_pdb.parent.mkdir(parents=True, exist_ok=True)
    try:
        unp_to_pdb = build_residue_map(uniprot, pdb_id, chain_id)
    except Exception as e:  # network/parse errors → fall through to copy
        log.warning("SIFTS lookup failed for %s_%s: %s — copying unchanged", pdb_id, chain_id, e)
        shutil.copy2(bronze_pdb, silver_pdb)
        return

    if not unp_to_pdb:
        log.warning("No SIFTS mapping for %s_%s — copying unchanged", pdb_id, chain_id)
        shutil.copy2(bronze_pdb, silver_pdb)
        return

    pdb_to_unp = {pdb_res: unp_pos for unp_pos, pdb_res in unp_to_pdb.items()}
    counts = renumber_pdb_chain(bronze_pdb, silver_pdb, chain_id, pdb_to_unp)
    log.info(
        "Renumbered %s chain %s: %d residue records shifted, %d unmapped, %d on other chains",
        pdb_id,
        chain_id,
        counts["renumbered"],
        counts["unmapped_in_chain"],
        counts["other"],
    )


def run_renumber(prefix: str, indir: Path, outdir: Path) -> None:
    """Process all structures listed in ``{prefix}_structures.json``."""
    metadata_in = indir / f"{prefix}_structures.json"
    if not metadata_in.exists():
        log.warning("%s missing — no structures to renumber.", metadata_in)
        outdir.mkdir(parents=True, exist_ok=True)
        return

    structures = json.loads(metadata_in.read_text())
    outdir.mkdir(parents=True, exist_ok=True)

    for s in structures:
        src = indir / s["file"]
        dst = outdir / s["file"]
        if not src.exists():
            log.warning("Bronze PDB missing: %s — skipping", src)
            continue

        if s["source"] == "alphafold":
            shutil.copy2(src, dst)
            log.info("AlphaFold %s already in UniProt numbering — copied", s["file"])
            continue

        renumber_one(
            bronze_pdb=src,
            silver_pdb=dst,
            uniprot=s["uniprot"],
            pdb_id=s["pdb_id"],
            chain_id=s["chain_id"],
        )

    # Mirror the metadata so analyze_site_models can find structures here
    shutil.copy2(metadata_in, outdir / metadata_in.name)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Step 2b: Renumber Bronze PDBs to UniProt numbering via SIFTS.",
    )
    add_common_args(ap)
    ap.add_argument(
        "--indir",
        type=Path,
        required=True,
        help="Bronze directory (output of fetch_structure.py).",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    ids = resolve_gene(args.gene_symbol, outdir=args.outdir, uniprot_override=args.uniprot)
    run_renumber(prefix=ids.prefix, indir=args.indir, outdir=ids.workdir)


if __name__ == "__main__":
    main()
