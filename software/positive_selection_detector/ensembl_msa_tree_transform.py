#!/usr/bin/env python3
"""Transform Bronze Ensembl downloads into a CodeML-ready alignment (Silver).

Reads the four artefacts produced by ``ensembl_msa_tree_extract.py`` from
``results/1_bronze/<prefix>/`` and produces a cleaned, length-filtered,
reference-species-restricted, TrimAl-cleaned CDS MSA, a pruned gene tree,
and a PHYLIP file ready for CodeML.

Pipeline (for prefix ``HLA_DQB1``):

    1. Copy/clean gene tree              → ``HLA_DQB1.tree``
    2. Filter AA/CDS for length parity   → ``HLA_DQB1_subset.aa.mafft.fasta``,
                                            ``HLA_DQB1_subset.cds.fasta``
    3. Filter to reference species panel → (same files, non-reference removed)
    4. Back-map CDS on the protein MSA
       via ``realign_nuc_on_aa.py``      → ``HLA_DQB1_subset.cds.mafft.fasta``
    5. Move target sequence to the top   → (same file, re-ordered)
    6. Clean with TrimAl                 → ``HLA_DQB1_subset.cds.mafft.trimal.fasta``,
                                            ``…trimal.html``, ``…trimal.cols``
    7. Prune gene tree with ``nw_prune`` → ``HLA_DQB1_subset.tree``
    8. Convert FASTA → PHYLIP via
       ``convert_fasta2phylip.py``       → ``HLA_DQB1_subset.cds.mafft.trimal.phy``

Outputs are written to ``results/2_silver/<prefix>/`` by default.

External tools required on ``PATH``:
    - ``trimal``   (https://github.com/inab/trimal)
    - ``nw_prune`` (Newick Utilities, https://gensoft.pasteur.fr/docs/newick-utils/)

The Compara REST endpoint already returns an aligned protein MSA, so MAFFT
is not required — the ``.aa.mafft.fasta`` and ``.cds.mafft.fasta`` filenames
are kept purely for compatibility with the tutorial's downstream commands.

Example:
    ensembl_msa_tree_transform.py --gene-symbol HLA-DQB1
    ensembl_msa_tree_transform.py --gene-symbol TRIM5 \\
                                  --indir results/1_bronze/TRIM5 \\
                                  --outdir results/2_silver/TRIM5
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path

from _common import (  # type: ignore
    add_common_args,
    check_tool,
    read_fasta_dict,
    resolve_gene,
    write_seq_fasta,
)
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ``reference_species`` (scripts/) is importable because importing
# ``_common`` puts that directory on sys.path.
from reference_species import is_reference_species  # type: ignore

# ── Repository layout ─────────────────────────────────────────────────
_REPO = Path(__file__).resolve().parent.parent.parent

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

_SCRIPTS = _REPO / "scripts"
_REALIGN_SCRIPT = _SCRIPTS / "realign_nuc_on_aa.py"
_CONVERT_SCRIPT = _SCRIPTS / "convert_fasta2phylip.py"

# Default layer directories (relative to the current working directory).
_DEFAULT_BRONZE_ROOT = Path("results/1_bronze")
_DEFAULT_SILVER_ROOT = Path("results/2_silver")


# =====================================================================
#  External-tool wrappers
# =====================================================================


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
    check_tool("trimal")
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
        "0.6",
        "-resoverlap",
        "0.65",
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
    check_tool("nw_prune")
    if not keep_ids:
        raise ValueError("Cannot prune tree: no sequences to keep.")
    log.info("Pruning tree to %d taxa", len(keep_ids))
    cmd = ["nw_prune", "-v", str(tree_file), *keep_ids]
    pruned = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.strip()
    with open(out_tree, "w") as fh:
        fh.write(f"{len(keep_ids)} 1\n")
        fh.write(pruned + "\n")


# =====================================================================
#  In-process helpers
# =====================================================================


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


def drop_empty_columns(aln: dict[str, str]) -> tuple[dict[str, str], int]:
    """Remove alignment columns that are all gaps. Returns (new_aln, n_dropped)."""
    if not aln:
        return aln, 0
    seqs = list(aln.values())
    width = len(seqs[0])
    keep = [i for i in range(width) if any(s[i] != "-" for s in seqs)]
    if len(keep) == width:
        return aln, 0
    trimmed = {sid: "".join(s[i] for i in keep) for sid, s in aln.items()}
    return trimmed, width - len(keep)


# =====================================================================
#  Pipeline
# =====================================================================


def run_transform(
    prefix: str,
    target_protein: str | None,
    bronze_dir: Path,
    silver_dir: Path,
) -> None:
    """Transform Bronze outputs into the CodeML-ready Silver outputs.

    Args:
        prefix: Filename stem, e.g. ``"HLA_DQB1"``.
        target_protein: Ensembl protein ID to place first in the MSA
            (the reference sequence). If ``None``, the MSA order is
            left unchanged.
        bronze_dir: Directory containing the four Bronze artefacts
            written by ``ensembl_msa_tree_extract.py``.
        silver_dir: Destination directory for all Silver outputs.
    """
    if not bronze_dir.is_dir():
        raise FileNotFoundError(f"Bronze input directory not found: {bronze_dir}")

    silver_dir.mkdir(parents=True, exist_ok=True)

    # Bronze inputs
    raw_tree = bronze_dir / f"{prefix}.nh"
    bronze_aa = bronze_dir / f"{prefix}.aa.fasta"
    bronze_cds = bronze_dir / f"{prefix}.cds.fasta"
    for path in (raw_tree, bronze_aa, bronze_cds):
        if not path.exists():
            raise FileNotFoundError(f"Required Bronze file missing: {path}")

    # Silver outputs — keep the tutorial filenames so downstream CodeML
    # commands run unchanged.
    clean_tree = silver_dir / f"{prefix}.tree"
    aa_msa = silver_dir / f"{prefix}_subset.aa.mafft.fasta"
    cds_subset = silver_dir / f"{prefix}_subset.cds.fasta"
    cds_msa = silver_dir / f"{prefix}_subset.cds.mafft.fasta"
    trimal_fasta = silver_dir / f"{prefix}_subset.cds.mafft.trimal.fasta"
    trimal_html = silver_dir / f"{prefix}_subset.cds.mafft.trimal.html"
    trimal_cols = silver_dir / f"{prefix}_subset.cds.mafft.trimal.cols"
    pruned_tree = silver_dir / f"{prefix}_subset.tree"
    phylip = silver_dir / f"{prefix}_subset.cds.mafft.trimal.phy"
    id_list = silver_dir / "id.list"

    # ── 1. Clean tree ────────────────────────────────────────────
    # The REST endpoint with ``nh_format=simple`` already returns bare
    # Ensembl IDs (no ``_species`` suffix), so no stripping is needed.
    # Running ``remove_ensembl_name_in_tree.py`` here would wrongly split
    # on every ``_`` and mangle native IDs such as ``MGP_CAROLIEiJ_P0044190``.
    clean_tree.write_text(raw_tree.read_text())

    # ── 2. Load Bronze AA + CDS, filter for length parity ────────
    aa_aligned = read_fasta_dict(bronze_aa)
    cds_map = read_fasta_dict(bronze_cds)

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

    # ── 3. Filter to reference species panel ─────────────────────
    # Drop non-reference species before back-mapping so that the
    # downstream realign/TrimAl steps only operate on the curated panel.
    before = len(consistent)
    consistent = {pid: pair for pid, pair in consistent.items() if is_reference_species(pid)}
    if not consistent:
        raise RuntimeError(
            "No sequences in reference species panel — check that "
            "--taxon includes the species families you need "
            "(e.g. 7742 = Vertebrata)."
        )
    log.info("Reference-species filter: kept %d / %d sequences", len(consistent), before)

    # Drop alignment columns that became all-gap after the species filter
    # so the back-mapped CDS MSA inherits a tighter, gap-free skeleton.
    aa_only = {pid: aln for pid, (aln, _) in consistent.items()}
    aa_only, dropped = drop_empty_columns(aa_only)
    if dropped:
        log.info("Dropped %d all-gap columns from the protein MSA", dropped)
    consistent = {pid: (aa_only[pid], cds) for pid, (_, cds) in consistent.items()}

    # Reorder so the target sequence leads every downstream alignment.
    if target_protein:
        if target_protein not in consistent:
            raise RuntimeError(
                f"Target protein {target_protein!r} was dropped by the "
                "length-parity or reference-species filter — cannot place "
                "it first in the alignment."
            )
        consistent = {
            target_protein: consistent[target_protein],
            **{pid: pair for pid, pair in consistent.items() if pid != target_protein},
        }
    else:
        log.warning("No --target protein supplied — leaving sequence order unchanged.")

    aa_records = [
        SeqRecord(Seq(aln), id=pid, description="") for pid, (aln, _) in consistent.items()
    ]
    SeqIO.write(aa_records, aa_msa, "fasta")
    write_seq_fasta({pid: cds for pid, (_, cds) in consistent.items()}, cds_subset)

    # ── 4. Back-map CDS onto the protein MSA ─────────────────────
    _run_python_script(
        _REALIGN_SCRIPT,
        str(aa_msa),
        str(cds_subset),
        str(cds_msa),
    )

    # ── 5. Re-assert target on top (safety net after realign) ────
    if target_protein:
        move_target_to_top(cds_msa, target_protein)

    # ── 6. Clean the MSA with TrimAl ─────────────────────────────
    run_trimal(cds_msa, trimal_fasta, trimal_html, trimal_cols)
    if not trimal_fasta.exists():
        log.warning(
            "TrimAl produced no output for %s (%s missing) — "
            "alignment likely too gappy or sequences too divergent. "
            "Skipping downstream Silver steps.",
            prefix,
            trimal_fasta,
        )
        sys.exit(0)
    if target_protein:
        move_target_to_top(trimal_fasta, target_protein)

    # ── 7. Prune the gene tree with nw_prune ─────────────────────
    keep_ids = ids_from_fasta(trimal_fasta)
    id_list.write_text("\n".join(keep_ids) + "\n")
    log.info("TrimAl kept %d sequences", len(keep_ids))
    run_nw_prune(clean_tree, keep_ids, pruned_tree)

    # ── 8. Convert cleaned FASTA → PHYLIP ────────────────────────
    _run_python_script(_CONVERT_SCRIPT, str(trimal_fasta), str(phylip))

    log.info("Done. Silver outputs in %s", silver_dir)
    log.info("  CDS MSA (trimmed) : %s", trimal_fasta)
    log.info("  CDS MSA (PHYLIP)  : %s", phylip)
    log.info("  Tree (pruned)     : %s", pruned_tree)


def run_transform_for_symbol(
    gene_symbol: str,
    outdir: Path | None,
    indir: Path | None,
    target: str | None = None,
    uniprot: str | None = None,
) -> None:
    """Resolve a gene symbol and run the Silver transform step."""
    ids = resolve_gene(gene_symbol, outdir=outdir, uniprot_override=uniprot)
    bronze_dir = indir if indir is not None else _DEFAULT_BRONZE_ROOT / ids.prefix
    silver_dir = outdir if outdir is not None else _DEFAULT_SILVER_ROOT / ids.prefix
    run_transform(
        prefix=ids.prefix,
        target_protein=target or ids.target,
        bronze_dir=bronze_dir,
        silver_dir=silver_dir,
    )


# =====================================================================
#  CLI
# =====================================================================


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    ap = argparse.ArgumentParser(
        description=(
            "Transform Bronze Ensembl downloads into a CodeML-ready "
            "alignment + tree (Silver layer at results/2_silver/<prefix>)."
        ),
    )
    add_common_args(ap)
    ap.add_argument(
        "--indir",
        type=Path,
        default=None,
        help="Bronze input directory. Defaults to results/1_bronze/<prefix>.",
    )
    ap.add_argument(
        "--target",
        default=None,
        help=(
            "Ensembl protein ID to place first in the MSA "
            "(default: the protein ID returned by the gene resolver)."
        ),
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    run_transform_for_symbol(
        gene_symbol=args.gene_symbol,
        outdir=args.outdir,
        indir=args.indir,
        target=args.target,
        uniprot=args.uniprot,
    )


if __name__ == "__main__":
    main()
