#!/usr/bin/env python3
"""Step 2-bis — Run HyPhy site models alongside CodeML M8.

HyPhy 2.5+ has no literal ``M8``; the closest analogues are:

    - **FUBAR**  — Bayesian per-site posterior of ω > 1
      (the closest analogue to CodeML M8's BEB table).
    - **BUSTED** — gene-wide test for episodic positive selection
      (parallel to the CodeML M8 vs M8a LRT).

This script runs FUBAR by default and optionally BUSTED.
It re-uses the alignment + tree produced by Step 1
(``prepare_data.py``); no extra preparation is needed.

Examples::

    run_hyphy_site_models.py --gene-symbol HLA-DQB1
    run_hyphy_site_models.py --gene-symbol HLA-DQB1 --method BUSTED
    run_hyphy_site_models.py --gene-symbol HLA-DQB1 --method BOTH

Outputs (under the gene workdir)::

    {prefix}_busted.json   — full BUSTED result (raw HyPhy JSON)
    {prefix}_busted.log    — captured stdout
    {prefix}_fubar.json    — full FUBAR result (raw HyPhy JSON)
    {prefix}_fubar.log     — captured stdout
    {prefix}_fubar_sites.tsv — positively-selected sites with
                               protein positions (when FUBAR is run)

External tools required on ``PATH``:
    - ``hyphy`` (>= 2.5)
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

from _common import (
    add_common_args,
    check_tool,
    log,
    map_sites_to_cds,
    resolve_gene,
)

METHODS = ("FUBAR", "BUSTED", "BOTH")
FUBAR_DEFAULT_THRESHOLD = 0.90


# =====================================================================
#  Locate Step 1 inputs
# =====================================================================


def locate_inputs(workdir: Path, prefix: str) -> tuple[Path, Path]:
    """Return ``(fasta, tree)`` produced by Step 1, or raise."""
    fasta = workdir / f"{prefix}_subset.cds.mafft.trimal.fasta"
    tree = workdir / f"{prefix}_subset.tree"
    if not fasta.exists() or not tree.exists():
        raise RuntimeError(
            f"Input files missing: {fasta}, {tree}. Run prepare_data.py (Step 1) first."
        )
    return fasta, tree


# =====================================================================
#  HyPhy invocation
# =====================================================================


def run_hyphy(
    method: str,
    fasta: Path,
    tree: Path,
    out_json: Path,
    log_path: Path,
) -> None:
    """Run a single HyPhy analysis (``busted`` or ``fubar``)."""
    if out_json.exists():
        log.info("Skipping HyPhy %s — %s already exists.", method, out_json.name)
        return

    check_tool("hyphy")
    cmd = [
        "hyphy",
        method.lower(),
        "--alignment",
        str(fasta),
        "--tree",
        str(tree),
        "--output",
        str(out_json),
    ]
    log.info("Running HyPhy %s (this can take some minutes)", method)
    with log_path.open("w") as fh:
        subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT)
    log.info("HyPhy %s finished → %s", method, out_json)


# =====================================================================
#  Result parsing / reporting
# =====================================================================


def report_busted(json_path: Path) -> None:
    """Print a one-line summary of the BUSTED result (LRT, p-value)."""
    if not json_path.exists():
        log.warning("BUSTED JSON missing: %s", json_path)
        return
    data = json.loads(json_path.read_text())
    test = data.get("test results", {})
    lrt = test.get("LRT")
    pval = test.get("p-value")

    print("\n=== HyPhy BUSTED ===")
    if lrt is None or pval is None:
        log.warning("Could not parse BUSTED test results from %s", json_path)
        return
    print(f"  LRT     = {lrt:10.4f}")
    print(f"  p-value = {pval:.3e}")


def _parse_fubar_table(json_path: Path) -> list[tuple[int, float, float, float]]:
    """Return ``[(trimal_pos, alpha, beta, prob_pos), …]`` from a FUBAR JSON."""
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


def write_fubar_sites_tsv(
    json_path: Path,
    tsv_path: Path,
    workdir: Path,
    prefix: str,
    target: str,
    uniprot: str | None,
    threshold: float = FUBAR_DEFAULT_THRESHOLD,
) -> None:
    """Write FUBAR positively-selected sites with target-protein positions.

    Mirrors the format of ``{prefix}_beb_sites.tsv``: only sites with
    ``Pr(beta > alpha) >= threshold`` are written.
    """
    rows = _parse_fubar_table(json_path)
    if not rows:
        return

    positives = [r for r in rows if r[3] >= threshold]
    full_aln = workdir / f"{prefix}_subset.cds.mafft.fasta"
    trimal_cols = workdir / f"{prefix}_subset.cds.mafft.trimal.cols"

    cds_mapping: dict[int, tuple[int | None, str]] = {}
    if positives and full_aln.exists() and trimal_cols.exists():
        cds_mapping = map_sites_to_cds(full_aln, trimal_cols, target, [r[0] for r in positives])
    elif positives:
        log.warning(
            "Cannot map trimal→CDS positions (missing %s or %s).",
            full_aln,
            trimal_cols,
        )

    lines: list[str] = []
    meta = []
    if target:
        meta.append(f"ensembl={target}")
    if uniprot:
        meta.append(f"uniprot={uniprot}")
    lines.append(
        "# FUBAR positively-selected sites (Pr(beta > alpha) >= "
        f"{threshold}). protein_pos is 1-based along the target protein "
        "sequence (" + "; ".join(meta) + ")."
    )
    lines.append("trimal_pos\tprotein_pos\taa\talpha\tbeta\tprob_pos")
    for site, alpha, beta, prob in positives:
        protein_pos, aa = cds_mapping.get(site, (None, ""))
        pos_str = "" if protein_pos is None else str(protein_pos)
        lines.append(f"{site}\t{pos_str}\t{aa}\t{alpha:.4f}\t{beta:.4f}\t{prob:.4f}")
    tsv_path.write_text("\n".join(lines) + "\n")
    log.info(
        "Wrote %d FUBAR sites (Pr >= %.2f) → %s",
        len(positives),
        threshold,
        tsv_path,
    )


def report_fubar(json_path: Path, threshold: float = FUBAR_DEFAULT_THRESHOLD) -> None:
    """Print a summary of FUBAR positively-selected sites above threshold."""
    if not json_path.exists():
        log.warning("FUBAR JSON missing: %s", json_path)
        return
    rows = _parse_fubar_table(json_path)
    if not rows:
        return
    n_pos = sum(1 for _, _, _, p in rows if p >= threshold)
    print("\n=== HyPhy FUBAR ===")
    print(f"  Sites with Pr(beta>alpha) >= {threshold}: {n_pos} / {len(rows)}")


# =====================================================================
#  Pipeline
# =====================================================================


def run_hyphy_pipeline(
    gene_symbol: str,
    method: str,
    outdir: Path | None,
    uniprot: str | None = None,
) -> None:
    """Run BUSTED and/or FUBAR on the Step-1 alignment/tree."""
    if method not in METHODS:
        raise ValueError(f"Unknown method {method!r}. Choose from: {METHODS}")

    ids = resolve_gene(gene_symbol, outdir=outdir, uniprot_override=uniprot)
    prefix = ids.prefix
    workdir = ids.workdir
    workdir.mkdir(parents=True, exist_ok=True)

    fasta, tree = locate_inputs(workdir, prefix)

    methods_to_run = ("FUBAR", "BUSTED") if method == "BOTH" else (method,)

    for m in methods_to_run:
        out_json = workdir / f"{prefix}_{m.lower()}.json"
        log_path = workdir / f"{prefix}_{m.lower()}.log"
        run_hyphy(m, fasta, tree, out_json, log_path)

    if "BUSTED" in methods_to_run:
        report_busted(workdir / f"{prefix}_busted.json")
    if "FUBAR" in methods_to_run:
        fubar_json = workdir / f"{prefix}_fubar.json"
        write_fubar_sites_tsv(
            fubar_json,
            workdir / f"{prefix}_fubar_sites.tsv",
            workdir=workdir,
            prefix=prefix,
            target=ids.target,
            uniprot=ids.uniprot,
        )
        report_fubar(fubar_json)


# =====================================================================
#  CLI
# =====================================================================


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Step 2-bis: Run HyPhy site models (BUSTED / FUBAR) "
        "alongside CodeML M8 on the Step-1 alignment.",
    )
    add_common_args(ap)
    ap.add_argument(
        "--method",
        choices=METHODS,
        default="FUBAR",
        help="Which HyPhy analysis to run (default: FUBAR — closest "
        "analogue to CodeML M8's BEB per-site posteriors).",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    run_hyphy_pipeline(
        gene_symbol=args.gene_symbol,
        method=args.method,
        outdir=args.outdir,
        uniprot=args.uniprot,
    )


if __name__ == "__main__":
    main()
