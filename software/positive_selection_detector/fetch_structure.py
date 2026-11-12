#!/usr/bin/env python3
"""Step 1b — Fetch 3D structures for the target protein.

Downloads:
    1. AlphaFold monomer for the UniProt accession.
    2. Top-N RCSB PDB structures covering the UniProt (SIFTS-ranked,
       greedy set-cover by residue span). Disable with ``--max-pdb 0``.
    3. Optionally a manually pinned PDB via ``--pdb`` (added on top).

Outputs (in ``--outdir``):
    - ``{prefix}_{uniprot}_AF.pdb``               – AlphaFold model
    - ``{prefix}_{PDB_ID}_{CHAIN}.pdb``           – one per RCSB structure
    - ``{prefix}_structures.json``                – list of structures with
                                                    metadata (source, pdb_id,
                                                    chain_id, uniprot, file)

Example:
    fetch_structure.py --gene-symbol HLA-DQB1
    fetch_structure.py --gene-symbol HLA-DQB1 --pdb 1UVQ
    fetch_structure.py --gene-symbol HLA-DQB1 --max-pdb 3
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

# ``sifts`` (lib/) is importable because importing ``_common`` puts that
# directory on sys.path.
from _common import (  # type: ignore
    ALPHAFOLD_API,
    ALPHAFOLD_FILES,
    RCSB_FILES,
    add_common_args,
    http_get,
    log,
    resolve_gene,
)
from sifts import get_sifts_mappings, pick_best_structures


def fetch_alphafold_pdb(uniprot: str, out_path: Path) -> bool:
    """Download the AlphaFold monomer PDB for a UniProt accession.

    Queries the AlphaFold prediction API first, then falls back to a
    v6 → v1 version sweep.
    """
    pdb_url: str | None = None
    api_raw = http_get(f"{ALPHAFOLD_API}/{uniprot}")
    if api_raw is not None:
        try:
            data = json.loads(api_raw)
        except json.JSONDecodeError:
            data = None
        if isinstance(data, list) and data:
            pdb_url = data[0].get("pdbUrl")

    if pdb_url:
        raw = http_get(pdb_url)
        if raw is not None:
            out_path.write_bytes(raw)
            log.info("AlphaFold model %s → %s", uniprot, out_path)
            return True

    for version in range(6, 0, -1):
        url = f"{ALPHAFOLD_FILES}/AF-{uniprot}-F1-model_v{version}.pdb"
        raw = http_get(url)
        if raw is not None:
            out_path.write_bytes(raw)
            log.info("AlphaFold model %s (v%d) → %s", uniprot, version, out_path)
            return True
    return False


def fetch_rcsb_pdb(pdb_id: str, out_path: Path) -> bool:
    """Download an experimental PDB structure from the RCSB."""
    url = f"{RCSB_FILES}/{pdb_id.upper()}.pdb"
    raw = http_get(url)
    if raw is None:
        return False
    out_path.write_bytes(raw)
    log.info("RCSB PDB %s → %s", pdb_id.upper(), out_path)
    return True


def _alphafold_filename(prefix: str, uniprot: str) -> str:
    return f"{prefix}_{uniprot}_AF.pdb"


def _rcsb_filename(prefix: str, pdb_id: str, chain_id: str) -> str:
    return f"{prefix}_{pdb_id.upper()}_{chain_id}.pdb"


def fetch_all_structures(
    outdir: Path,
    prefix: str,
    uniprot: str | None,
    pinned_pdb: str | None,
    max_pdb: int,
) -> list[dict]:
    """Fetch AlphaFold + top-N RCSB PDBs (SIFTS) for ``uniprot``.

    Returns a list of structure metadata dicts:
        {source, pdb_id, chain_id, uniprot, file, unp_start, unp_end}
    where ``file`` is the basename of the saved PDB.
    """
    structures: list[dict] = []

    if not uniprot:
        log.warning(
            "Could not resolve a UniProt accession for %s — skipping all structure downloads.",
            prefix,
        )
        return structures

    # 1. AlphaFold
    af_name = _alphafold_filename(prefix, uniprot)
    if fetch_alphafold_pdb(uniprot, outdir / af_name):
        structures.append(
            {
                "source": "alphafold",
                "pdb_id": "AF",
                "chain_id": "A",
                "uniprot": uniprot,
                "file": af_name,
                "unp_start": None,
                "unp_end": None,
            }
        )

    # 2. SIFTS-ranked RCSB PDBs
    sifts_pdbs: list[dict] = []
    if max_pdb > 0:
        log.info("Querying SIFTS for PDB structures of %s …", uniprot)
        mappings = get_sifts_mappings(uniprot)
        sifts_pdbs = pick_best_structures(mappings, max_n=max_pdb)
        log.info(
            "SIFTS picked %d structures: %s",
            len(sifts_pdbs),
            ", ".join(f"{m['pdb_id']}_{m['chain_id']}" for m in sifts_pdbs),
        )

    # 3. Honour --pdb override (added if not already in SIFTS picks)
    selected = list(sifts_pdbs)
    if pinned_pdb:
        pid_lc = pinned_pdb.lower()
        if not any(m["pdb_id"].lower() == pid_lc for m in selected):
            # Look up the chain via the full SIFTS list (all_mappings); fall back
            # to chain "A" if not found.
            full = get_sifts_mappings(uniprot)
            chain = next(
                (m["chain_id"] for m in full if m["pdb_id"].lower() == pid_lc),
                "A",
            )
            selected.append(
                {
                    "pdb_id": pinned_pdb,
                    "chain_id": chain,
                    "unp_start": None,
                    "unp_end": None,
                }
            )

    for m in selected:
        pdb_id = m["pdb_id"].upper()
        chain = m["chain_id"]
        fname = _rcsb_filename(prefix, pdb_id, chain)
        if fetch_rcsb_pdb(pdb_id, outdir / fname):
            structures.append(
                {
                    "source": "rcsb",
                    "pdb_id": pdb_id,
                    "chain_id": chain,
                    "uniprot": uniprot,
                    "file": fname,
                    "unp_start": m.get("unp_start"),
                    "unp_end": m.get("unp_end"),
                }
            )

    return structures


def run_fetch_structure(
    gene_symbol: str,
    outdir: Path | None,
    uniprot: str | None = None,
    pdb_id: str | None = None,
    max_pdb: int = 5,
) -> None:
    """Fetch AlphaFold + RCSB PDBs for a gene and write a structures.json."""
    ids = resolve_gene(gene_symbol, outdir=outdir, uniprot_override=uniprot)
    ids.workdir.mkdir(parents=True, exist_ok=True)

    structures = fetch_all_structures(
        outdir=ids.workdir,
        prefix=ids.prefix,
        uniprot=ids.uniprot,
        pinned_pdb=pdb_id,
        max_pdb=max_pdb,
    )

    metadata_path = ids.workdir / f"{ids.prefix}_structures.json"
    metadata_path.write_text(json.dumps(structures, indent=2) + "\n")
    log.info("Wrote %s (%d structures)", metadata_path, len(structures))

    if not structures:
        log.warning("No 3D structures could be fetched.")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Step 1b: Fetch 3D structures (AlphaFold + SIFTS-ranked RCSB PDBs).",
    )
    add_common_args(ap)
    ap.add_argument(
        "--pdb",
        default=None,
        help="Pin this RCSB PDB accession in addition to SIFTS picks.",
    )
    ap.add_argument(
        "--max-pdb",
        type=int,
        default=5,
        help="Maximum number of RCSB PDBs from SIFTS (0 to disable, default: 5).",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    run_fetch_structure(
        gene_symbol=args.gene_symbol,
        outdir=args.outdir,
        uniprot=args.uniprot,
        pdb_id=args.pdb,
        max_pdb=args.max_pdb,
    )


if __name__ == "__main__":
    main()
