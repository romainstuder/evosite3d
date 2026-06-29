#!/usr/bin/env python3
"""Step 3 — Analyse site-model results.

Reads the per-site selection outputs produced by Step 2
(``run_codeml_site_models.py``) and/or Step 2-bis
(``run_hyphy_site_models.py``) and generates visualisation files:

    - ``{prefix}_beb.png``       — per-site dN/dS Manhattan plot (CodeML)
    - ``{prefix}_beb.jlv``       — Jalview annotation (CodeML BEB + dN/dS)
    - ``{prefix}_fubar.jlv``     — Jalview annotation (HyPhy FUBAR)
    - ``{prefix}_*.pml``         — PyMOL script highlighting selected sites

CodeML BEB and HyPhy FUBAR are handled independently: whichever inputs
are present drive the corresponding outputs. A HyPhy-only run (no CodeML)
still produces ``{prefix}_fubar.jlv`` and a FUBAR-only PyMOL script.

Assumes Step 1 (``prepare_data.py``) has fetched the 3D structure, and
that Step 2 produced ``{prefix}_beb.txt`` / ``{prefix}_beb_sites.tsv``
and/or Step 2-bis produced ``{prefix}_fubar.json`` /
``{prefix}_fubar_sites.tsv``.

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

from _common import (  # type: ignore
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


def _load_trimal_aa_columns(trimal_file_path: Path) -> list[int]:
    """Return the unique 1-based AA column indices kept by TrimAl.

    Reads the ``#ColumnsMap`` line from a TrimAl ``.cols`` file (CDS column
    indices), converts each to its amino-acid position, and de-duplicates
    while preserving order.
    """
    column_str_list: list[str] = []
    with open(trimal_file_path) as file_in:
        for line in file_in:
            line = line.rstrip()
            if line.startswith("#"):
                cleaned_line = line.replace(",", "")
                column_str_list = cleaned_line.split()[1:]

    column_int_list = [int(x) for x in column_str_list]
    trimal_aa_cols = [(pos // 3) + 1 for pos in column_int_list]
    return list(dict.fromkeys(trimal_aa_cols))


def _track_over_columns(trimal_cols: list[int], value_by_trimpos: dict[int, float]) -> list[float]:
    """Spread per-trimmed-column values over the full alignment columns.

    ``trimal_cols`` is the ordered list of original alignment AA columns
    kept by TrimAl; its position in that list (1-based) is the trimmed
    column index used as the key in ``value_by_trimpos``. Returns a value
    for every alignment column up to the last kept one (0 elsewhere).
    """
    corr_site_dict = {aln_col: i + 1 for i, aln_col in enumerate(trimal_cols)}
    out: list[float] = []
    for site in range(1, trimal_cols[-1] + 1):
        trim_pos = corr_site_dict.get(site)
        out.append(value_by_trimpos.get(trim_pos, 0.0) if trim_pos is not None else 0.0)
    return out


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

    trimal_cols = _load_trimal_aa_columns(trimal_file_path)

    beb_probs_dict: dict[int, float] = {}
    mean_omegas_dict: dict[int, float] = {}
    for row in rows:
        cols = row.split()
        aln_pos = int(cols[0])
        beb_probs_dict[aln_pos] = float(cols[12])
        mean_omegas_dict[aln_pos] = float(cols[14])

    beb_probs = _track_over_columns(trimal_cols, beb_probs_dict)
    mean_omegas = _track_over_columns(trimal_cols, mean_omegas_dict)

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


def _load_fubar_table(json_path: Path) -> list[tuple[int, float, float, float]]:
    """Return ``[(trimal_pos, alpha, beta, prob_pos), …]`` from a FUBAR JSON.

    ``trimal_pos`` is the 1-based codon position in the trimmed alignment,
    matching the ``trimal_pos`` column of ``{prefix}_fubar_sites.tsv``.
    """
    data = json.loads(json_path.read_text())
    mle = data.get("MLE", {})
    headers = [h[0] for h in mle.get("headers", [])]
    content = mle.get("content", {}).get("0", [])
    if not headers or not content or "Prob[alpha<beta]" not in headers:
        log.warning("Could not parse FUBAR result table from %s", json_path)
        return []
    i_a = headers.index("alpha")
    i_b = headers.index("beta")
    i_p = headers.index("Prob[alpha<beta]")
    return [(site, row[i_a], row[i_b], row[i_p]) for site, row in enumerate(content, 1)]


def write_fubar_jalview_annotation(
    path: Path,
    trimal_file_path: Path,
    fubar_json: Path,
    threshold: float = 0.90,
) -> None:
    """Write a Jalview annotation file from HyPhy FUBAR results.

    Mirrors :func:`write_jalview_annotation` but for FUBAR, producing:
        - **FUBAR Pr(beta>alpha)** — bar graph of posterior probability,
          teal above ``threshold``.
        - **dN/dS (beta/alpha)** — line graph with ω = 1 reference
          (capped at 99 to keep the line graph finite when alpha ≈ 0).
    """
    rows = _load_fubar_table(fubar_json)
    if not rows:
        log.warning("Empty FUBAR table — skipping FUBAR Jalview annotation.")
        return

    trimal_cols = _load_trimal_aa_columns(trimal_file_path)

    prob_dict: dict[int, float] = {}
    omega_dict: dict[int, float] = {}
    for trimal_pos, alpha, beta, prob in rows:
        prob_dict[trimal_pos] = prob
        omega_dict[trimal_pos] = min(beta / max(alpha, 1e-3), 99.0)

    probs = _track_over_columns(trimal_cols, prob_dict)
    omegas = _track_over_columns(trimal_cols, omega_dict)

    prob_vals: list[str] = []
    for p in probs:
        colour = "00808a" if p > threshold else "000000"  # teal above threshold
        prob_vals.append(f"{p:.4f},{p:.4f},{colour}")

    omega_vals = [f"{w:.4f}" for w in omegas]

    lines = [
        "JALVIEW_ANNOTATION",
        f"BAR_GRAPH\tFUBAR Pr(beta>alpha)\t"
        f"FUBAR posterior probability of positive selection (beta > alpha)"
        f"\t{'|'.join(prob_vals)}",
        f"LINE_GRAPH\tdN/dS (beta/alpha)\t"
        f"FUBAR posterior beta/alpha per site"
        f"\t{'|'.join(omega_vals)}",
        "GRAPHLINE\tdN/dS (beta/alpha)\t1.0\tneutral\t000000",
    ]
    path.write_text("\n".join(lines) + "\n")
    log.info("Wrote FUBAR Jalview annotation (%d sites) → %s", len(probs), path)


# =====================================================================
#  PyMOL script
# =====================================================================


def _parse_sites_tsv(tsv_path: Path, prob_col: int) -> list[tuple[int | None, float]]:
    """Read a sites TSV and return ``(protein_pos, prob)`` pairs.

    ``{prefix}_beb_sites.tsv`` and ``{prefix}_fubar_sites.tsv`` share the same
    ``trimal_pos``/``protein_pos``/… column layout; ``prob_col`` selects the
    probability column (3 for CodeML BEB, 5 for HyPhy FUBAR).
    """
    results: list[tuple[int | None, float]] = []
    for line in tsv_path.read_text().splitlines():
        if line.startswith("#") or line.startswith("trimal_pos"):
            continue
        cols = line.split("\t")
        protein_pos = int(cols[1]) if cols[1] else None
        results.append((protein_pos, float(cols[prob_col])))
    return results


def write_pymol_script(
    pml: Path,
    pdb_path: Path,
    sites_tsv: Path | None,
    chain: str = "A",
    resi_offset: int = 0,
    threshold_high: float = 0.95,
    threshold_low: float = 0.50,
    fubar_sites_tsv: Path | None = None,
    fubar_threshold: float = 0.90,
) -> None:
    """Write a PyMOL ``.pml`` that highlights positively-selected sites.

    CodeML BEB layers (when ``sites_tsv`` is given):
        - yellow spheres for ``Pr(w > 1) >= threshold_high``
        - orange sticks for ``threshold_low <= Pr(w > 1) < threshold_high``

    If a HyPhy FUBAR sites file is provided, an extra layer of cyan
    sticks is drawn first for FUBAR sites with
    ``Pr(beta > alpha) >= fubar_threshold``. CodeML layers are drawn
    after, so overlapping sites display the higher-confidence colour.

    ``sites_tsv`` may be ``None`` (HyPhy-only run), in which case only the
    FUBAR layer is drawn.
    """
    high, low = [], []
    if sites_tsv is not None:
        for protein_pos, prob in _parse_sites_tsv(sites_tsv, prob_col=3):
            if protein_pos is None:
                continue
            if prob >= threshold_high:
                high.append(protein_pos)
            elif prob >= threshold_low:
                low.append(protein_pos)

    fubar: list[int] = []
    if fubar_sites_tsv is not None and fubar_sites_tsv.exists():
        for protein_pos, prob in _parse_sites_tsv(fubar_sites_tsv, prob_col=5):
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
    fubar_json = workdir / f"{prefix}_fubar.json"
    fubar_sites_tsv = workdir / f"{prefix}_fubar_sites.tsv"
    trimal_cols = workdir / f"{prefix}_subset.cds.mafft.trimal.cols"

    # CodeML BEB: dN/dS plot + Jalview annotation (require beb.txt from Step 2)
    if beb_txt.exists():
        beb_png = workdir / f"{prefix}_beb.png"
        jlv = workdir / f"{prefix}_beb.jlv"
        try:
            plot_beb_per_site(beb_txt, beb_png, threshold=plot_threshold)
            write_jalview_annotation(jlv, trimal_cols, beb_txt, threshold=plot_threshold)
        except (RuntimeError, subprocess.CalledProcessError) as e:
            log.warning("Per-site dN/dS plot failed: %s", e)
    else:
        log.warning(
            "%s missing — skipping dN/dS plot and CodeML Jalview annotation. "
            "Run run_codeml_site_models.py (Step 2) first.",
            beb_txt,
        )

    # HyPhy FUBAR: Jalview annotation (requires fubar.json from Step 2-bis)
    if fubar_json.exists() and trimal_cols.exists():
        write_fubar_jalview_annotation(workdir / f"{prefix}_fubar.jlv", trimal_cols, fubar_json)

    # PyMOL script: needs a structure plus at least one of BEB / FUBAR sites.
    beb_arg = beb_sites_tsv if beb_sites_tsv.exists() else None
    fubar_arg = fubar_sites_tsv if fubar_sites_tsv.exists() else None
    if beb_arg is None and fubar_arg is None:
        log.warning(
            "Neither %s nor %s present — skipping PyMOL script. Run "
            "run_codeml_site_models.py (Step 2) or run_hyphy_site_models.py "
            "(Step 2-bis) first.",
            beb_sites_tsv,
            fubar_sites_tsv,
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

    taxon_tag = f"_{taxon}" if taxon is not None else ""

    for s in structures:
        pml_name = f"{prefix}{taxon_tag}_{s['pdb_id']}.pml"
        # The Silver step renumbered each chain to UniProt → no offset needed
        # unless the caller explicitly overrides via --resi-offset.
        struct_chain = s.get("chain_id") or chain
        write_pymol_script(
            workdir / pml_name,
            s["path"],
            beb_arg,
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
