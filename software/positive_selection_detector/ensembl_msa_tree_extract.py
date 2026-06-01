#!/usr/bin/env python3
"""Fetch Ensembl Compara gene tree, protein MSA, and CDS sequences.

Downloads the raw artefacts that the downstream
``ensembl_msa_tree_transform.py`` script consumes:

    1. Compara gene tree (Newick)        → ``<prefix>.nh``
    2. Leaf-name list                    → ``<prefix>_names.txt``
    3. Compara protein MSA (aligned)     → ``<prefix>.aa.fasta``
    4. Per-protein CDS sequences         → ``<prefix>.cds.fasta``

Outputs are written to ``results/1_bronze/<prefix>/`` by default.

Example:
    ensembl_msa_tree_extract.py --gene-symbol HLA-DQB1
    ensembl_msa_tree_extract.py --gene-symbol TRIM5 --taxon 7742
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

# ``ensembl`` (lib/) and ``reference_species`` (scripts/) are importable
# because importing ``_common`` puts those directories on sys.path.
from _common import (  # type: ignore
    add_common_args,
    resolve_gene,
    write_seq_fasta,
)
from ensembl import (  # type: ignore
    extract_protein_ids_from_newick,
    fetch_cds_sequences,
    fetch_compara_protein_msa,
    fetch_gene_tree_newick,
)
from reference_species import is_reference_species  # type: ignore

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

# Default Bronze-layer output root (relative to the current working directory).
_DEFAULT_ROOT = Path("results/1_bronze")


# =====================================================================
#  Helpers
# =====================================================================


def write_leaf_names(newick: str, out_path: Path) -> list[str]:
    """Extract Ensembl protein IDs from a Newick tree and write one per line."""
    ids = extract_protein_ids_from_newick(newick)
    out_path.write_text("\n".join(ids) + "\n")
    return ids


# =====================================================================
#  Pipeline
# =====================================================================


def run_extract(
    ensembl_gene: str,
    prefix: str,
    outdir: Path,
    taxon: int = 9347,
) -> None:
    """Download raw Ensembl tree + protein MSA + CDS into ``outdir``.

    Args:
        ensembl_gene: Ensembl gene stable ID (e.g. ``"ENSG00000179344"``).
        prefix: Filename stem, e.g. ``"HLA_DQB1"``.
        outdir: Destination directory for the four outputs.
        taxon: NCBI taxon ID used to prune the Compara gene tree
            (default 9347 = Boreoeutheria).
    """
    outdir.mkdir(parents=True, exist_ok=True)

    raw_tree = outdir / f"{prefix}.nh"
    names_txt = outdir / f"{prefix}_names.txt"
    aa_fasta = outdir / f"{prefix}.aa.fasta"
    cds_fasta = outdir / f"{prefix}.cds.fasta"

    if raw_tree.exists() and aa_fasta.exists() and cds_fasta.exists():
        log.info("Found existing outputs in %s – skipping fetch.", outdir)
        return

    # ── 1. Fetch gene tree ───────────────────────────────────────
    log.info("Fetching Compara gene tree for %s (taxon=%d)", ensembl_gene, taxon)
    newick = fetch_gene_tree_newick(ensembl_gene, prune_taxon=taxon)
    if not newick:
        raise RuntimeError(f"Failed to download gene tree for {ensembl_gene}")
    raw_tree.write_text(newick + "\n")

    # ── 2. Leaf names (mirrors the tutorial's `nw_labels -I` step) ──
    leaf_ids = write_leaf_names(newick, names_txt)
    log.info("Gene tree contains %d leaves", len(leaf_ids))

    # ── 3. Fetch protein MSA from Compara (keyed by protein ID) ──
    log.info("Fetching Compara protein MSA")
    aa_aligned = fetch_compara_protein_msa(ensembl_gene, prune_taxon=taxon)
    if not aa_aligned:
        raise RuntimeError(f"Failed to download protein MSA for {ensembl_gene}")
    aa_aligned = {pid: aa_aligned[pid] for pid in aa_aligned if is_reference_species(pid)}
    log.info("Protein MSA has %d sequences", len(aa_aligned))
    write_seq_fasta(aa_aligned, aa_fasta)

    # ── 4. Fetch CDS for every protein in the MSA ────────────────
    aa_ids = list(aa_aligned.keys())
    cds_map = fetch_cds_sequences(ensembl_gene, aa_ids)
    if not cds_map:
        raise RuntimeError("No CDS sequences could be fetched from Ensembl.")
    write_seq_fasta(cds_map, cds_fasta)

    log.info("Done. Outputs in %s", outdir)
    log.info("  Gene tree (raw) : %s", raw_tree)
    log.info("  Protein MSA     : %s", aa_fasta)
    log.info("  CDS sequences   : %s", cds_fasta)


def run_extract_for_symbol(
    gene_symbol: str,
    outdir: Path | None,
    taxon: int,
    uniprot: str | None = None,
) -> None:
    """Resolve a gene symbol and run the extract step."""
    ids = resolve_gene(gene_symbol, outdir=outdir, uniprot_override=uniprot)
    target_dir = outdir if outdir is not None else _DEFAULT_ROOT / ids.prefix
    run_extract(
        ensembl_gene=ids.gene,
        prefix=ids.prefix,
        outdir=target_dir,
        taxon=taxon,
    )


# =====================================================================
#  CLI
# =====================================================================


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    ap = argparse.ArgumentParser(
        description=(
            "Fetch raw Ensembl Compara gene tree, protein MSA, and CDS "
            "sequences into the layer (i.e. results/1_bronze/<prefix>)."
        ),
    )
    add_common_args(ap)
    ap.add_argument(
        "--taxon",
        type=int,
        default=9347,
        help="NCBI taxon ID to prune the Compara gene tree "
        "(default 9347 = Boreoeutheria; 9443 = Primates, "
        "40674 = Mammalia, 7742 = Vertebrata).",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    run_extract_for_symbol(
        gene_symbol=args.gene_symbol,
        outdir=args.outdir,
        taxon=args.taxon,
        uniprot=args.uniprot,
    )


if __name__ == "__main__":
    main()
