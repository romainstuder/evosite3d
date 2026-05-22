#!/usr/bin/env python3
"""Step 3 — Analyse CodeML site-model results.

Reads the BEB outputs produced by Step 2 (``run_codeml_site_models.py``)
and generates visualisation files:

    - ``{prefix}_beb.png``       — per-site dN/dS Manhattan plot
    - ``{prefix}_beb.jlv``       — Jalview annotation (BEB + dN/dS tracks)
    - ``{prefix}_sites.pml``     — PyMOL script highlighting BEB sites

Assumes that Step 1 (``prepare_data.py``) has fetched the 3D structure
and Step 2 has produced ``{prefix}_beb.txt`` and
``{prefix}_beb_sites.tsv``.

Example:
    python analyse_site_models.py --gene-symbol HLA-DQB1
    python analyse_site_models.py --gene-symbol HLA-DQB1 --pdb 1UVQ --chain B --resi-offset 32
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from ._common import (
    PLOT_SCRIPT,
    add_common_args,
    log,
    resolve_gene,
)

# =====================================================================
#  Per-site dN/dS plot
# =====================================================================


def plot_beb_per_site(beb_txt: Path, out_png: Path, threshold: float = 0.50) -> None:
    """Run ``plot_dnds_per_site.py`` to produce the dN/dS plot."""
    cmd = [
        sys.executable,
        str(PLOT_SCRIPT),
        str(beb_txt),
        "-o",
        str(out_png),
        "-t",
        str(threshold),
        "--no-show",
    ]
    log.info("Plotting dN/dS per site → %s", out_png)
    subprocess.run(cmd, check=True)


# =====================================================================
#  Jalview annotation
# =====================================================================


def write_jalview_annotation(
    path: Path,
    trimal_file_path: Path,
    beb_txt: Path,
    threshold: float = 0.50,
) -> None:
    """Write a Jalview annotation file with BEB and dN/dS tracks.

    Produces two annotation tracks:
        - **BEB Pr(w>1)** — bar graph of posterior probability.
        - **dN/dS (posterior mean)** — line graph with ω = 1 reference.
    """
    rows = beb_txt.read_text().splitlines()
    if not rows:
        log.warning("Empty BEB table — skipping Jalview annotation.")
        return

    # Load first sequences
    # first_aa_seq = next(SeqIO.parse(msa_aa_aln, format="fasta"))
    # print("seq_id", first_aa_seq.id)
    # print("seq", first_aa_seq.seq)

    # Load trimal cols
    def load_trimal_columns(trimal_file_path):
        column_str_list = []
        with open(trimal_file_path, "r") as file_in:
            for line in file_in:
                line = line.rstrip()
                if line.startswith("#"):
                    cleaned_line = line.replace(",", "")
                    column_str_list = cleaned_line.split()[1:]

        column_int_list = [int(x) for x in column_str_list]
        trimal_aa_cols = [(pos // 3) + 1 for pos in column_int_list]
        return list(dict.fromkeys(trimal_aa_cols))

    trimal_cols = load_trimal_columns(trimal_file_path)

    beb_probs_dict: dict[int, float] = {}
    mean_omegas_dict: dict[int, float] = {}
    for row in rows:
        cols = row.split()
        aln_pos = int(cols[0])
        beb_probs_dict[aln_pos] = float(cols[12])
        mean_omegas_dict[aln_pos] = float(cols[14])

    beb_probs: list[float] = []
    mean_omegas: list[float] = []
    corr_site_dict = {aln_pos: i + 1 for i, aln_pos in enumerate(trimal_cols)}

    for i in range(trimal_cols[-1]):
        site = i + 1
        if site in corr_site_dict:
            aln_pos = corr_site_dict[site]
            if aln_pos in beb_probs_dict:
                beb_probs.append(beb_probs_dict[aln_pos])
                mean_omegas.append(mean_omegas_dict[aln_pos])
            else:
                beb_probs.append(0)
                mean_omegas.append(0)
        else:
            beb_probs.append(0)
            mean_omegas.append(0)

    beb_vals: list[str] = []
    for p in beb_probs:
        colour = "ff0000" if p > threshold else "000000"
        beb_vals.append(f"{p:.4f},{p:.4f},{colour}")

    omega_vals: list[str] = []
    for w in mean_omegas:
        omega_vals.append(f"{w:.4f}")

    lines = [
        "JALVIEW_ANNOTATION",
        f"BAR_GRAPH\tBEB Pr(w>1)\t"
        f"BEB posterior probability of positive selection (M8 site model)"
        f"\t{'|'.join(beb_vals)}",
        f"LINE_GRAPH\tdN/dS (posterior mean)\t"
        f"Posterior mean omega per site (M8 BEB)"
        f"\t{'|'.join(omega_vals)}",
        "GRAPHLINE\tdN/dS (posterior mean)\t1.0\tneutral\t000000",
    ]
    path.write_text("\n".join(lines) + "\n")
    log.info("Wrote Jalview annotation (%d sites) → %s", len(beb_probs), path)


# =====================================================================
#  PyMOL script
# =====================================================================


def _parse_beb_sites_tsv(tsv_path: Path) -> list[tuple[int | None, float]]:
    """Read ``{prefix}_beb_sites.tsv`` and return (protein_pos, prob) pairs."""
    results: list[tuple[int | None, float]] = []
    for line in tsv_path.read_text().splitlines():
        if line.startswith("#") or line.startswith("trimal_pos"):
            continue
        cols = line.split("\t")
        protein_pos = int(cols[1]) if cols[1] else None
        prob = float(cols[3])
        results.append((protein_pos, prob))
    return results


def _parse_fubar_sites_tsv(tsv_path: Path) -> list[tuple[int | None, float]]:
    """Read ``{prefix}_fubar_sites.tsv`` and return (protein_pos, prob) pairs."""
    results: list[tuple[int | None, float]] = []
    for line in tsv_path.read_text().splitlines():
        if line.startswith("#") or line.startswith("trimal_pos"):
            continue
        cols = line.split("\t")
        protein_pos = int(cols[1]) if cols[1] else None
        prob = float(cols[5])
        results.append((protein_pos, prob))
    return results


def write_pymol_script(
    pml: Path,
    pdb_path: Path,
    sites_tsv: Path,
    chain: str = "A",
    resi_offset: int = 0,
    threshold_high: float = 0.95,
    threshold_low: float = 0.50,
    fubar_sites_tsv: Path | None = None,
    fubar_threshold: float = 0.90,
) -> None:
    """Write a PyMOL ``.pml`` that highlights BEB-positive sites.

    CodeML BEB layers:
        - yellow spheres for ``Pr(w > 1) >= threshold_high``
        - orange sticks for ``threshold_low <= Pr(w > 1) < threshold_high``

    If a HyPhy FUBAR sites file is provided, an extra layer of cyan
    sticks is drawn first for FUBAR sites with
    ``Pr(beta > alpha) >= fubar_threshold``. CodeML layers are drawn
    after, so overlapping sites display the higher-confidence colour.
    """
    parsed = _parse_beb_sites_tsv(sites_tsv)

    high, low = [], []
    for protein_pos, prob in parsed:
        if protein_pos is None:
            continue
        if prob >= threshold_high:
            high.append(protein_pos)
        elif prob >= threshold_low:
            low.append(protein_pos)

    fubar: list[int] = []
    if fubar_sites_tsv is not None and fubar_sites_tsv.exists():
        for protein_pos, prob in _parse_fubar_sites_tsv(fubar_sites_tsv):
            if protein_pos is None or prob < fubar_threshold:
                continue
            fubar.append(protein_pos)

    high_expr = "+".join(str(r) for r in high) if high else ""
    low_expr = "+".join(str(r) for r in low) if low else ""
    fubar_expr = "+".join(str(r) for r in fubar) if fubar else ""

    pdb_abs = pdb_path.resolve()
    lines = [
        "# PyMOL visualisation of positively-selected sites.",
        "# CodeML M8 BEB sites: yellow spheres (>=0.95), orange sticks (>=0.50).",
        "# HyPhy FUBAR sites:   cyan sticks (Pr(beta>alpha) >= 0.90).",
        "# Residues numbered by position in the target protein sequence",
        "# (Ensembl translation / UniProt canonical isoform).",
        f"# Generated by {Path(__file__).name}.",
        f"load {pdb_abs}, structure",
        f"select target, structure and chain {chain}",
        "hide everything",
        "show cartoon, structure",
        "colour grey70, structure",
    ]
    if resi_offset:
        lines.append(
            f"alter target, resi=str(int(resi)+{resi_offset})  "
            f"# align crystal numbering to target protein sequence"
        )
        lines.append("rebuild")
    if fubar_expr:
        lines += [
            f"select sites_FUBAR{int(fubar_threshold * 100)}, target and resi {fubar_expr}",
            f"show sticks, sites_FUBAR{int(fubar_threshold * 100)}",
            f"colour cyan, sites_FUBAR{int(fubar_threshold * 100)}",
        ]
    if low_expr:
        lines += [
            f"select sites_BEB{int(threshold_low * 100)}, target and resi {low_expr}",
            f"show sticks, sites_BEB{int(threshold_low * 100)}",
            f"colour orange, sites_BEB{int(threshold_low * 100)}",
        ]
    if high_expr:
        lines += [
            f"select sites_BEB{int(threshold_high * 100)}, target and resi {high_expr}",
            f"show spheres, sites_BEB{int(threshold_high * 100)}",
            f"colour yellow, sites_BEB{int(threshold_high * 100)}",
        ]
    lines += [
        "bg_color white",
        "util.cnc",
        "zoom target",
    ]
    pml.write_text("\n".join(lines) + "\n")
    log.info(
        "Wrote PyMOL script → %s (%d BEB sphere, %d BEB stick, %d FUBAR stick)",
        pml,
        len(high),
        len(low),
        len(fubar),
    )


# =====================================================================
#  Find the structure file in the workdir
# =====================================================================


def find_structures(workdir: Path, prefix: str) -> list[dict]:
    """Return all structures available in ``workdir``.

    Reads ``{prefix}_structures.json`` if present (written by
    ``fetch_structure.py``). Each item is the metadata dict augmented
    with ``path`` (the PDB on disk).

    Falls back to globbing ``*.pdb`` for backwards compatibility.
    """
    metadata = workdir / f"{prefix}_structures.json"
    if metadata.exists():
        items = json.loads(metadata.read_text())
        out = []
        for s in items:
            path = workdir / s["file"]
            if path.exists():
                out.append({**s, "path": path})
        return out

    # Legacy fallback — single AlphaFold PDB only
    candidates = sorted(workdir.glob(f"{prefix}_*_AF.pdb"))
    if not candidates:
        candidates = sorted(workdir.glob("*.pdb"))
    return [
        {
            "source": "unknown",
            "pdb_id": p.stem.split("_")[-1] or "PDB",
            "chain_id": "A",
            "uniprot": None,
            "file": p.name,
            "path": p,
        }
        for p in candidates
    ]


# =====================================================================
#  Pipeline
# =====================================================================


def run_analysis(
    gene_symbol: str,
    outdir: Path | None,
    uniprot: str | None = None,
    pdb_id: str | None = None,
    chain: str = "A",
    resi_offset: int = 0,
    plot_threshold: float = 0.50,
    taxon: int | None = None,
) -> None:
    """Generate Jalview annotation, dN/dS plot and PyMOL script."""
    ids = resolve_gene(gene_symbol, outdir=outdir, uniprot_override=uniprot)
    prefix = ids.prefix
    workdir = ids.workdir

    if not workdir.exists():
        raise RuntimeError(
            f"Working directory {workdir} does not exist. "
            "Run prepare_data.py (Step 1) and run_codeml_site_models.py (Step 2) first."
        )

    beb_txt = workdir / f"{prefix}_beb.txt"
    beb_sites_tsv = workdir / f"{prefix}_beb_sites.tsv"

    # dN/dS plot + Jalview annotation (require beb.txt from Step 2)
    if beb_txt.exists():
        beb_png = workdir / f"{prefix}_beb.png"
        jlv = workdir / f"{prefix}_beb.jlv"
        # msa_aa_aln= workdir / f"{prefix}_subset.aa.mafft.fasta"
        trimal_cols = workdir / f"{prefix}_subset.cds.mafft.trimal.cols"
        try:
            plot_beb_per_site(beb_txt, beb_png, threshold=plot_threshold)
            write_jalview_annotation(jlv, trimal_cols, beb_txt, threshold=plot_threshold)
        except (RuntimeError, subprocess.CalledProcessError) as e:
            log.warning("Per-site dN/dS plot failed: %s", e)
    else:
        log.warning(
            "%s missing — skipping dN/dS plot and Jalview annotation. "
            "Run run_codeml_site_models.py (Step 2) first.",
            beb_txt,
        )

    # PyMOL script (requires beb_sites.tsv from Step 2 + structure from Step 1)
    if not beb_sites_tsv.exists():
        log.warning(
            "%s missing — skipping PyMOL script. Run run_codeml_site_models.py (Step 2) first.",
            beb_sites_tsv,
        )
        return

    structures = find_structures(workdir, prefix)
    if pdb_id:
        pid_lc = pdb_id.lower()
        structures = [s for s in structures if s["pdb_id"].lower() == pid_lc]
    if not structures:
        log.warning(
            "No structure file found in %s — skipping PyMOL script. "
            "Run fetch_structure.py (Step 1b) or pass --pdb.",
            workdir,
        )
        return

    fubar_sites_tsv = workdir / f"{prefix}_fubar_sites.tsv"
    fubar_arg = fubar_sites_tsv if fubar_sites_tsv.exists() else None
    taxon_tag = f"_{taxon}" if taxon is not None else ""

    for s in structures:
        pml_name = f"{prefix}{taxon_tag}_{s['pdb_id']}.pml"
        # The Silver step renumbered each chain to UniProt → no offset needed
        # unless the caller explicitly overrides via --resi-offset.
        struct_chain = s.get("chain_id") or chain
        write_pymol_script(
            workdir / pml_name,
            s["path"],
            beb_sites_tsv,
            chain=struct_chain,
            resi_offset=resi_offset,
            fubar_sites_tsv=fubar_arg,
        )
        log.info("PyMOL script for %s_%s → %s", s["pdb_id"], struct_chain, workdir / pml_name)


# =====================================================================
#  CLI
# =====================================================================


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Step 3: Analyse CodeML results — dN/dS plot, "
        "Jalview annotation, PyMOL visualisation.",
    )
    add_common_args(ap)
    ap.add_argument(
        "--pdb",
        default=None,
        help="PDB accession to locate in the workdir for the PyMOL script.",
    )
    ap.add_argument(
        "--chain",
        default="A",
        help="Chain ID of the target in the structure (default A).",
    )
    ap.add_argument(
        "--resi-offset",
        type=int,
        default=0,
        help="Residue-numbering offset for the PyMOL script "
        "(e.g. 32 for the tutorial's 1UVQ example).",
    )
    ap.add_argument(
        "--plot-threshold",
        type=float,
        default=0.50,
        help="BEB probability threshold for the dN/dS plot (default 0.50).",
    )
    ap.add_argument(
        "--taxon",
        type=int,
        default=None,
        help="NCBI taxon ID — embedded in PML filenames as {prefix}_{taxon}_{pdb_id}.pml.",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    run_analysis(
        gene_symbol=args.gene_symbol,
        outdir=args.outdir,
        uniprot=args.uniprot,
        pdb_id=args.pdb,
        chain=args.chain,
        resi_offset=args.resi_offset,
        plot_threshold=args.plot_threshold,
        taxon=args.taxon,
    )


if __name__ == "__main__":
    main()
