#!/usr/bin/env python3
"""Fetch Ensembl CDS orthologues + Compara gene tree and prepare CodeML input.

Automates the manual download and preparation steps of the tutorial
``tutorials/pervasive_selection_site_models`` and preserves the same
filename nomenclature so the downstream CodeML commands from the
tutorial run unchanged.

Pipeline (for gene ``ENSG00000179344`` with prefix ``HLA_DQB1``):

    1. Fetch Compara gene tree          → ``HLA_DQB1.nh``
    2. Strip ``_ENS…`` species tags     → ``HLA_DQB1.tree``
    3. Fetch protein MSA from Compara   → ``HLA_DQB1_subset.aa.mafft.fasta``
    4. Fetch per-protein CDS sequences  → ``HLA_DQB1.cds.fasta``,
                                           ``HLA_DQB1_subset.cds.fasta``
    5. Back-map CDS on the protein MSA
       via ``realign_nuc_on_aa.py``     → ``HLA_DQB1_subset.cds.mafft.fasta``
    6. Move target sequence to the top  → (same file, re-ordered)
    6b. Filter to reference species     → (same file, non-reference species removed)
    7. Clean with TrimAl                → ``HLA_DQB1_subset.cds.mafft.trimal.fasta``
                                           ``…trimal.html``, ``…trimal.cols``
    8. Prune gene tree with ``nw_prune`` → ``HLA_DQB1_subset.tree``
    9. Convert FASTA → PHYLIP via
       ``convert_fasta2phylip.py``      → ``HLA_DQB1_subset.cds.mafft.trimal.phy``

External tools required on ``PATH``:
    - ``trimal``    (https://github.com/inab/trimal)
    - ``nw_prune`` (Newick Utilities, https://gensoft.pasteur.fr/docs/newick-utils/)

The Compara REST endpoint already returns an aligned protein MSA, so
MAFFT is not required — we keep the ``.aa.mafft.fasta`` filename purely
for compatibility with the tutorial's downstream commands.

Example:
    fetch_ensembl_msa_tree.py --gene ENSG00000179344 \\
                              --prefix HLA_DQB1 \\
                              --target ENSP00000407332 \\
                              --outdir HLA_DQB1_workdir
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Make the protein_characteriser and positive_selection_detector modules importable
_REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO / "lib"))
sys.path.insert(0, str(_REPO / "software" / "positive_selection_detector"))

from ensembl import (  # noqa: E402
    extract_protein_ids_from_newick,
    fetch_cds_sequences,
    fetch_compara_protein_msa,
    fetch_gene_tree_newick,
)
from reference_species import is_reference_species  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


_SCRIPTS = _REPO / "scripts"
_REALIGN_SCRIPT = _SCRIPTS / "realign_nuc_on_aa.py"
_CONVERT_SCRIPT = _SCRIPTS / "convert_fasta2phylip.py"


# =====================================================================
#  External-tool wrappers
# =====================================================================


def _check_tool(name: str) -> None:
    """Raise ``RuntimeError`` if an external binary is missing from ``PATH``."""
    if shutil.which(name) is None:
        raise RuntimeError(f"Required tool '{name}' not found on PATH.")


def _run_python_script(script: Path, *args: str, stdout: Path | None = None) -> None:
    """Run one of the repository's helper Python scripts via subprocess."""
    cmd = [sys.executable, str(script), *args]
    log.info("$ %s", " ".join(cmd))
    if stdout is not None:
        with open(stdout, "w") as fh:
            subprocess.run(cmd, check=True, stdout=fh)
    else:
        subprocess.run(cmd, check=True)


def run_trimal(in_fasta: Path, out_fasta: Path, html: Path, cols: Path) -> None:
    """Clean an MSA with TrimAl using the tutorial's settings."""
    _check_tool("trimal")
    log.info("Running TrimAl on %s", in_fasta)
    cmd = [
        "trimal",
        "-in",
        str(in_fasta),
        "-out",
        str(out_fasta),
        "-htmlout",
        str(html),
        "-gt",
        "0.7",
        "-resoverlap",
        "0.75",
        "-seqoverlap",
        "85",
        "-colnumbering",
        "-block",
        "3",
    ]
    with open(cols, "w") as colfh:
        subprocess.run(cmd, check=True, stdout=colfh)


def run_nw_prune(tree_file: Path, keep_ids: list[str], out_tree: Path) -> None:
    """Prune a Newick tree to ``keep_ids`` and prepend the PAML header line."""
    _check_tool("nw_prune")
    if not keep_ids:
        raise ValueError("Cannot prune tree: no sequences to keep.")
    log.info("Pruning tree to %d taxa", len(keep_ids))
    cmd = ["nw_prune", "-v", str(tree_file), *keep_ids]
    pruned = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.strip()
    with open(out_tree, "w") as fh:
        fh.write(f"{len(keep_ids)} 1\n")
        fh.write(pruned + "\n")


# =====================================================================
#  In-process helpers (small enough to keep inline)
# =====================================================================


def write_cds_fasta(cds: dict[str, str], out_fasta: Path) -> None:
    """Write a ``{id: sequence}`` mapping as multi-FASTA."""
    records = [SeqRecord(Seq(seq), id=pid, description="") for pid, seq in cds.items()]
    SeqIO.write(records, out_fasta, "fasta")


def subset_cds_fasta(in_fasta: Path, keep_ids: list[str], out_fasta: Path) -> None:
    """Write only the records whose ID is in ``keep_ids`` (preserves order)."""
    wanted = set(keep_ids)
    records = [r for r in SeqIO.parse(in_fasta, "fasta") if r.id in wanted]
    SeqIO.write(records, out_fasta, "fasta")


def write_leaf_names(newick: str, out_path: Path) -> list[str]:
    """Extract Ensembl protein IDs from a Newick tree and write one per line."""
    ids = extract_protein_ids_from_newick(newick)
    out_path.write_text("\n".join(ids) + "\n")
    return ids


def move_target_to_top(fasta: Path, target: str) -> None:
    """Re-order a FASTA file so the record matching ``target`` comes first.

    CodeML uses the first sequence as the positional reference, so the
    human/target protein must lead the alignment.
    """
    records = list(SeqIO.parse(fasta, "fasta"))
    target_rec = next((r for r in records if r.id == target), None)
    if target_rec is None:
        log.warning("Target %s not found in %s — leaving order unchanged", target, fasta)
        return
    reordered = [target_rec] + [r for r in records if r.id != target]
    SeqIO.write(reordered, fasta, "fasta")
    log.info("Moved %s to top of %s", target, fasta)


def ids_from_fasta(fasta: Path) -> list[str]:
    """Return the list of record IDs (in file order) from a FASTA file."""
    return [r.id for r in SeqIO.parse(fasta, "fasta")]


# =====================================================================
#  Pipeline
# =====================================================================


def run_pipeline(
    ensembl_gene: str,
    prefix: str,
    target_protein: str | None,
    outdir: Path,
    taxon: int = 9347,
) -> None:
    """Run the full CDS-level preparation pipeline for CodeML.

    Args:
        ensembl_gene: Ensembl gene stable ID (e.g. ``"ENSG00000179344"``).
        prefix: Filename stem, e.g. ``"HLA_DQB1"``.
        target_protein: Ensembl protein ID to place first in the MSA
            (the reference sequence). If ``None``, the human ID returned
            by Compara is used.
        outdir: Directory for all output files.
        taxon: NCBI taxon ID used to prune the Compara gene tree
            (default 9347 = Boreoeutheria, matching the tutorial's scope).
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # Tutorial-nomenclature filenames
    raw_tree = outdir / f"{prefix}.nh"
    clean_tree = outdir / f"{prefix}.tree"
    names_txt = outdir / f"{prefix}_names.txt"
    aa_msa = outdir / f"{prefix}_subset.aa.mafft.fasta"
    cds_fasta = outdir / f"{prefix}.cds.fasta"
    cds_subset = outdir / f"{prefix}_subset.cds.fasta"
    cds_msa = outdir / f"{prefix}_subset.cds.mafft.fasta"
    trimal_fasta = outdir / f"{prefix}_subset.cds.mafft.trimal.fasta"
    trimal_html = outdir / f"{prefix}_subset.cds.mafft.trimal.html"
    trimal_cols = outdir / f"{prefix}_subset.cds.mafft.trimal.cols"
    pruned_tree = outdir / f"{prefix}_subset.tree"
    phylip = outdir / f"{prefix}_subset.cds.mafft.trimal.phy"
    id_list = outdir / "id.list"

    # ── 1. Fetch gene tree ───────────────────────────────────────
    log.info("Fetching Compara gene tree for %s (taxon=%d)", ensembl_gene, taxon)
    newick = fetch_gene_tree_newick(ensembl_gene, prune_taxon=taxon)
    if not newick:
        raise RuntimeError(f"Failed to download gene tree for {ensembl_gene}")
    raw_tree.write_text(newick + "\n")

    # ── 2. Clean tree ────────────────────────────────────────────
    # The REST endpoint with ``nh_format=simple`` already returns bare
    # Ensembl IDs (no ``_species`` suffix), so no stripping is needed.
    # Running ``remove_ensembl_name_in_tree.py`` here would wrongly split
    # on every ``_`` and mangle native IDs such as ``MGP_CAROLIEiJ_P0044190``.
    clean_tree.write_text(newick + "\n")

    # Leaf names (mirrors the tutorial's `nw_labels -I` step)
    leaf_ids = write_leaf_names(newick, names_txt)
    log.info("Gene tree contains %d leaves", len(leaf_ids))

    # ── 3. Fetch protein MSA from Compara (keyed by protein ID) ──
    log.info("Fetching Compara protein MSA")
    aa_aligned = fetch_compara_protein_msa(ensembl_gene, prune_taxon=taxon)
    if not aa_aligned:
        raise RuntimeError(f"Failed to download protein MSA for {ensembl_gene}")
    log.info("Protein MSA has %d sequences", len(aa_aligned))

    # ── 4. Fetch CDS for every protein in the MSA ────────────────
    aa_ids = list(aa_aligned.keys())
    cds_map = fetch_cds_sequences(ensembl_gene, aa_ids)
    if not cds_map:
        raise RuntimeError("No CDS sequences could be fetched from Ensembl.")

    # Drop sequences whose CDS length is incompatible with the protein
    # sequence (len_cds in {3·len_prot, 3·len_prot + 3}). realign_nuc_on_aa
    # errors out on any mismatch, so filter them here to keep the pipeline
    # running on the well-behaved majority of orthologues.
    consistent: dict[str, tuple[str, str]] = {}
    for pid, aligned in aa_aligned.items():
        cds = cds_map.get(pid)
        if not cds:
            continue
        aa_len = len(aligned.replace("-", ""))
        if len(cds) not in (aa_len * 3, aa_len * 3 + 3):
            log.debug("Skipping %s: CDS %d nt vs protein %d aa", pid, len(cds), aa_len)
            continue
        consistent[pid] = (aligned, cds)
    log.info(
        "Kept %d / %d orthologues with compatible AA + CDS lengths",
        len(consistent),
        len(aa_aligned),
    )
    if not consistent:
        raise RuntimeError("No orthologues with matching AA / CDS lengths.")

    # Write the AA MSA (realign_nuc_on_aa will read it back)
    aa_records = [
        SeqRecord(Seq(aln), id=pid, description="") for pid, (aln, _) in consistent.items()
    ]
    SeqIO.write(aa_records, aa_msa, "fasta")

    # Write full + subset CDS FASTAs (tutorial uses both filenames)
    write_cds_fasta({pid: cds for pid, cds in cds_map.items()}, cds_fasta)
    write_cds_fasta({pid: cds for pid, (_, cds) in consistent.items()}, cds_subset)

    # ── 5. Back-map CDS onto the protein MSA ─────────────────────
    _run_python_script(
        _REALIGN_SCRIPT,
        str(aa_msa),
        str(cds_subset),
        str(cds_msa),
    )

    # ── 6. Move the target/human sequence to the top ─────────────
    if target_protein:
        move_target_to_top(cds_msa, target_protein)
    else:
        log.warning("No --target protein supplied — leaving sequence order unchanged.")

    # ── 6b. Filter to reference species panel ────────────────────
    # Remove non-reference species *before* TrimAl so that gap and
    # overlap statistics reflect only the curated panel.
    all_records = list(SeqIO.parse(cds_msa, "fasta"))
    ref_records = [r for r in all_records if is_reference_species(r.id)]
    if not ref_records:
        raise RuntimeError(
            "No sequences in reference species panel — check that "
            "--taxon includes the species families you need "
            "(e.g. 7742 = Vertebrata)."
        )
    log.info(
        "Reference-species filter: kept %d / %d sequences",
        len(ref_records),
        len(all_records),
    )
    SeqIO.write(ref_records, cds_msa, "fasta")

    # ── 7. Clean the MSA with TrimAl ─────────────────────────────
    run_trimal(cds_msa, trimal_fasta, trimal_html, trimal_cols)

    # ── 8. Prune the gene tree with nw_prune ─────────────────────
    keep_ids = ids_from_fasta(trimal_fasta)
    id_list.write_text("\n".join(keep_ids) + "\n")
    log.info("TrimAl kept %d sequences", len(keep_ids))
    run_nw_prune(clean_tree, keep_ids, pruned_tree)

    # ── 9. Convert cleaned FASTA → PHYLIP ────────────────────────
    _run_python_script(_CONVERT_SCRIPT, str(trimal_fasta), str(phylip))

    log.info("Done. Outputs in %s", outdir)
    log.info("  CDS MSA (trimmed) : %s", trimal_fasta)
    log.info("  CDS MSA (PHYLIP)  : %s", phylip)
    log.info("  Tree (pruned)     : %s", pruned_tree)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    ap = argparse.ArgumentParser(
        description="Fetch Ensembl CDS MSA + gene tree and prepare them for CodeML.",
    )
    ap.add_argument(
        "--gene",
        required=True,
        help="Ensembl gene stable ID (e.g. ENSG00000179344).",
    )
    ap.add_argument(
        "--prefix",
        required=True,
        help="Filename stem used for all outputs (e.g. HLA_DQB1).",
    )
    ap.add_argument(
        "--target",
        default=None,
        help="Ensembl protein ID to place first in the MSA "
        "(default: the human sequence returned by Compara).",
    )
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path("."),
        help="Output directory.",
    )
    ap.add_argument(
        "--taxon",
        type=int,
        default=9347,
        help="NCBI taxon ID to prune the Compara gene tree "
        "(default 9347 = Boreoeutheria, matching the tutorial; "
        "9443 = Primates, 40674 = Mammalia, 7742 = Vertebrata).",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    run_pipeline(args.gene, args.prefix, args.target, args.outdir, taxon=args.taxon)


if __name__ == "__main__":
    main()
