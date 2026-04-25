#!/usr/bin/env python3
"""Step 2 — Run CodeML site models and extract results.

Supports two modes of operation:

**Run a single model** (designed for Nextflow parallelism)::

    run_codeml_site_models.py --gene-symbol HLA-DQB1 --model M8
    run_codeml_site_models.py --gene-symbol HLA-DQB1 --model M8a

**Extract BEB results** (after both models have finished)::

    run_codeml_site_models.py --gene-symbol HLA-DQB1 --extract-only

**Run everything sequentially** (standalone, no Nextflow)::

    run_codeml_site_models.py --gene-symbol HLA-DQB1

Assumes that Step 1 (``prepare_data.py``) has already produced the
PHYLIP alignment and pruned tree in the working directory.

Outputs produced:
    - ``{prefix}_M8.ctl``, ``{prefix}_M8a.ctl``  — control files
    - ``{prefix}_M8.mlc``, ``{prefix}_M8a.mlc``   — CodeML results
    - ``{prefix}_M8.rst.txt``                      — preserved rst
    - ``{prefix}_beb_sites.tsv``  — BEB sites with protein positions
    - ``{prefix}_beb.txt``        — cleaned per-site BEB probability table

External tools required on ``PATH``:
    - ``codeml`` (PAML)
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from _common import (
    GET_POS_SCRIPT,
    BebSite,
    LnL,
    add_common_args,
    check_tool,
    log,
    parse_beb_sites,
    parse_lnl,
    resolve_gene,
)
from scipy.stats import chi2

# =====================================================================
#  Model definitions
# =====================================================================


# Each model is defined by its NSsites value and fix_omega setting.
@dataclass
class ModelDef:
    """CodeML site-model definition."""

    nssites: str
    fix_omega: int


MODELS: dict[str, ModelDef] = {
    "M8": ModelDef(nssites="8", fix_omega=0),
    "M8a": ModelDef(nssites="8", fix_omega=1),
}


# =====================================================================
#  Control-file templates
# =====================================================================


_CTL_COMMON = """\
     seqfile = {seqfile}  * sequence data file name
    treefile = {treefile}  * tree structure file name
     outfile = {outfile}  * main result file name

       noisy = 9   * 0,1,2,3,9: how much rubbish on the screen
     verbose = 1   * 1: detailed output, 0: concise output
     runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

     seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
   CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
       clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
      aaDist = 0   * 0:equal, +:geometric; -:linear, {{1-5:G1974,Miyata,c,p,v}}
       model = 0   * models for codons:
                   * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
     NSsites = {nssites}   * 0:one w; 1:NearlyNeutral; 2:PositiveSelection;
                   * 3:discrete; 4:freqs; 5:gamma;6:2gamma;
                   * 7:beta;8:beta&w;9:beta&gamma;10:3normal
       icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
       Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all
   fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
       kappa = 2   * initial or fixed kappa
   fix_omega = {fix_omega}   * 1: omega or omega_1 fixed, 0: estimate
       omega = 1   * initial or fixed omega, for codons or codon-based AAs
       getSE = 0       * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0       * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
  Small_Diff = .45e-6  * Default value.
   cleandata = 0       * remove sites with ambiguity data (1:yes, 0:no)?
 fix_blength = 0       * 0: ignore, -1: random, 1: initial, 2: fixed
"""


def write_ctl(
    ctl_path: Path,
    seqfile: str,
    treefile: str,
    outfile: str,
    nssites: str,
    fix_omega: int,
) -> None:
    """Write a CodeML control file."""
    ctl_path.write_text(
        _CTL_COMMON.format(
            seqfile=seqfile,
            treefile=treefile,
            outfile=outfile,
            nssites=nssites,
            fix_omega=fix_omega,
        )
    )
    log.info("Wrote %s", ctl_path)


# =====================================================================
#  Running CodeML
# =====================================================================


def run_codeml(ctl: Path, workdir: Path, mlc: Path | None = None) -> None:
    """Run ``codeml`` on a control file in an isolated subdirectory.

    Creates a subdirectory named after the CTL stem, symlinks inputs,
    runs ``codeml``, and moves outputs back.
    """
    if mlc is not None and mlc.exists():
        log.info("Skipping codeml on %s – %s already exists.", ctl.name, mlc.name)
        return
    check_tool("codeml")

    rundir = workdir / ctl.stem
    rundir.mkdir(exist_ok=True)

    for src in workdir.iterdir():
        if src.suffix in (".phy", ".tree") and src.is_file():
            dst = rundir / src.name
            if not dst.exists():
                dst.symlink_to(src.resolve())

    shutil.copy2(ctl, rundir / ctl.name)

    log.info("Running codeml on %s (this can take some minutes)", ctl.name)
    subprocess.run(["codeml", ctl.name], check=True, cwd=rundir)

    mlc_names = {f for f in rundir.iterdir() if f.suffix == ".mlc"}
    for path in sorted(rundir.iterdir()):
        if path.is_symlink() or path.name == ctl.name:
            continue
        if path in mlc_names:
            dest = workdir / path.name
        else:
            dest = workdir / f"{ctl.stem}_{path.name}"
        shutil.move(str(path), str(dest))

    for path in rundir.iterdir():
        path.unlink()
    rundir.rmdir()


# =====================================================================
#  LRT reporting
# =====================================================================


def lrt(name: str, null: LnL, alt: LnL) -> str:
    """Format a single likelihood-ratio test line."""
    stat = 2.0 * (alt.lnl - null.lnl)
    df = alt.np - null.np
    p = chi2.sf(stat, df) if df > 0 and stat > 0 else float("nan")
    return f"{name:>7}  2ΔlnL = {stat:10.4f}   df = {df}   p = {p:.3e}"


def report_lrts(prefix: str, outdir: Path) -> None:
    """Parse CodeML outputs and print the M8a-M8 LRT."""
    m8_mlc = outdir / f"{prefix}_M8.mlc"
    m8a_mlc = outdir / f"{prefix}_M8a.mlc"
    if not m8_mlc.exists() or not m8a_mlc.exists():
        log.warning("Missing mlc files – cannot compute LRTs (%s, %s)", m8_mlc, m8a_mlc)
        return

    m8_entries = parse_lnl(m8_mlc)
    m8a_entries = parse_lnl(m8a_mlc)
    if not m8_entries or not m8a_entries:
        log.warning(
            "No lnL entries found: %d in M8, %d in M8a",
            len(m8_entries),
            len(m8a_entries),
        )
        return

    m8 = m8_entries[-1]
    m8a = m8a_entries[-1]

    print("\n=== Log-likelihoods ===")
    for name, e in [("M8", m8), ("M8a", m8a)]:
        print(f"  {name:>3}  np = {e.np:3d}   lnL = {e.lnl:12.4f}")

    print("\n=== Likelihood-ratio test ===")
    print(lrt("M8a-M8", m8a, m8))
    print()


# =====================================================================
#  BEB site mapping and extraction
# =====================================================================


_GET_POS_RE = re.compile(r"site_index:\s*(\d+)\s+CDS pos:\s*(\d+|N/A)\s+AA:\s*(\S+)")


def map_sites_to_cds(
    full_aln: Path,
    trimal_cols: Path,
    target: str,
    sites: list[int],
) -> dict[int, tuple[int | None, str]]:
    """Map CodeML/trimal column indices to CDS positions in the target."""
    if not sites:
        return {}
    site_arg = " ".join(str(s) for s in sites)
    cmd = [
        sys.executable,
        str(GET_POS_SCRIPT),
        str(full_aln),
        str(trimal_cols),
        target,
        site_arg,
    ]
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    mapping: dict[int, tuple[int | None, str]] = {}
    for line in result.stdout.splitlines():
        m = _GET_POS_RE.search(line)
        if not m:
            continue
        pos = int(m.group(1))
        cds = None if m.group(2) == "N/A" else int(m.group(2))
        mapping[pos] = (cds, m.group(3))
    return mapping


def write_beb_sites_tsv(
    path: Path,
    sites: list[BebSite],
    cds_mapping: dict[int, tuple[int | None, str]],
    target: str | None = None,
    uniprot: str | None = None,
) -> None:
    """Write a tab-separated summary of BEB sites + target-protein positions."""
    lines: list[str] = []
    if target or uniprot:
        meta = []
        if target:
            meta.append(f"ensembl={target}")
        if uniprot:
            meta.append(f"uniprot={uniprot}")
        lines.append(
            "# BEB-positive sites (M8 site model). protein_pos is 1-based "
            "along the target protein sequence (" + "; ".join(meta) + ")."
        )
    lines.append("trimal_pos\tprotein_pos\taa\tprob\tsignif\tmean_w\tse")
    for s in sites:
        protein_pos, _ = cds_mapping.get(s.trimal_pos, (None, s.aa))
        pos_str = "" if protein_pos is None else str(protein_pos)
        lines.append(
            f"{s.trimal_pos}\t{pos_str}\t{s.aa}\t{s.prob:.3f}\t"
            f"{s.stars}\t{s.mean_w:.3f}\t{s.se:.3f}"
        )
    path.write_text("\n".join(lines) + "\n")
    log.info("Wrote %d BEB sites → %s", len(sites), path)


# =====================================================================
#  BEB probability table extraction from rst
# =====================================================================


_BEB_TABLE_HEADER_RE = re.compile(
    r"(Bayes Empirical Bayes \(BEB\)|Naive Empirical Bayes \(NEB\))"
    r" probabilities for (\d+) classes"
)
_BEB_TABLE_ROW_RE = re.compile(r"^\s*\d+\s+\S\s+")


def extract_beb_table(rst_txt: Path, out_beb_txt: Path) -> int:
    """Extract the per-site posterior-probabilities block from ``rst``.

    Locates the BEB (or NEB fallback) table, normalises spacing, and
    pads to 11 classes so the positive-selection class sits at column 12.

    Returns:
        The number of per-site rows written.
    """
    text = rst_txt.read_text()

    beb_match = None
    neb_match = None
    for m in _BEB_TABLE_HEADER_RE.finditer(text):
        if m.group(1).startswith("Bayes"):
            beb_match = m
        else:
            neb_match = m

    match = beb_match or neb_match
    if match is None:
        raise RuntimeError(
            "No 'Bayes Empirical Bayes (BEB) probabilities for N classes' "
            f"or NEB fallback block found in {rst_txt} — cannot build plot."
        )

    n_classes = int(match.group(2))
    which = "BEB" if beb_match is match else "NEB"

    rows: list[str] = []
    seen_data = False
    for line in text[match.start() :].splitlines()[1:]:
        if _BEB_TABLE_ROW_RE.match(line):
            seen_data = True
            rows.append(line.replace("( ", "("))
        elif seen_data and not line.strip():
            break

    if not rows:
        raise RuntimeError(f"{which} probability table was empty in {rst_txt}")

    if n_classes < 11:
        pad_cols = 11 - n_classes
        pad = "0.0 " * pad_cols
        padded: list[str] = []
        for row in rows:
            parts = row.split(None, 2)
            if len(parts) < 3:
                padded.append(row)
                continue
            padded.append(f"{parts[0]} {parts[1]} {pad}{parts[2]}")
        rows = padded

    out_beb_txt.write_text("\n".join(rows) + "\n")
    log.info(
        "Wrote %d-row %s table (%d classes%s) → %s",
        len(rows),
        which,
        n_classes,
        ", padded to 11 for plot layout" if n_classes < 11 else "",
        out_beb_txt,
    )
    return len(rows)


# =====================================================================
#  Single-model runner (for Nextflow)
# =====================================================================


def run_single_model(
    gene_symbol: str,
    model_name: str,
    outdir: Path | None,
    uniprot: str | None = None,
) -> None:
    """Write the CTL file for one model and run codeml."""
    if model_name not in MODELS:
        raise ValueError(f"Unknown model {model_name!r}. Choose from: {list(MODELS)}")

    ids = resolve_gene(gene_symbol, outdir=outdir, uniprot_override=uniprot)
    ids.workdir.mkdir(parents=True, exist_ok=True)

    prefix = ids.prefix
    workdir = ids.workdir
    model = MODELS[model_name]

    phy = workdir / f"{prefix}_subset.cds.mafft.trimal.phy"
    tree = workdir / f"{prefix}_subset.tree"

    if not phy.exists() or not tree.exists():
        raise RuntimeError(
            f"Input files missing: {phy}, {tree}. Run prepare_data.py (Step 1) first."
        )

    ctl = workdir / f"{prefix}_{model_name}.ctl"
    mlc_name = f"{prefix}_{model_name}.mlc"

    write_ctl(
        ctl,
        seqfile=phy.name,
        treefile=tree.name,
        outfile=mlc_name,
        nssites=model.nssites,
        fix_omega=model.fix_omega,
    )

    mlc_path = workdir / mlc_name
    run_codeml(ctl, workdir, mlc=mlc_path)

    # Preserve rst for M8 (needed for BEB extraction)
    if model_name == "M8":
        rst_src = workdir / f"{prefix}_M8_rst"
        rst_dst = workdir / f"{prefix}_M8.rst.txt"
        if rst_src.exists() and (
            not rst_dst.exists() or rst_src.stat().st_mtime >= rst_dst.stat().st_mtime
        ):
            shutil.move(str(rst_src), str(rst_dst))
            log.info("Moved %s → %s", rst_src.name, rst_dst.name)


# =====================================================================
#  Extract results (after both models have run)
# =====================================================================


def run_extract(
    gene_symbol: str,
    outdir: Path | None,
    uniprot: str | None = None,
) -> None:
    """Report LRT, extract BEB sites and BEB probability table.

    Requires that both M8 and M8a ``.mlc`` files already exist in the
    working directory.
    """
    ids = resolve_gene(gene_symbol, outdir=outdir, uniprot_override=uniprot)
    prefix = ids.prefix
    workdir = ids.workdir

    report_lrts(prefix, workdir)

    # Extract BEB sites and map to CDS/protein positions
    mlc = workdir / f"{prefix}_M8.mlc"
    rst = workdir / f"{prefix}_M8.rst.txt"
    full_aln = workdir / f"{prefix}_subset.cds.mafft.fasta"
    trimal_cols = workdir / f"{prefix}_subset.cds.mafft.trimal.cols"

    if not mlc.exists():
        log.warning("%s missing — cannot extract BEB sites.", mlc)
        return

    sites = parse_beb_sites(mlc)
    if not sites:
        log.warning("No BEB sites parsed from %s.", mlc.name)
        return
    log.info("Parsed %d BEB sites from %s", len(sites), mlc.name)

    cds_mapping: dict[int, tuple[int | None, str]] = {}
    if full_aln.exists() and trimal_cols.exists():
        cds_mapping = map_sites_to_cds(
            full_aln,
            trimal_cols,
            ids.target,
            [s.trimal_pos for s in sites],
        )
    else:
        log.warning(
            "Cannot map trimal→CDS positions (missing %s or %s).",
            full_aln,
            trimal_cols,
        )

    write_beb_sites_tsv(
        workdir / f"{prefix}_beb_sites.tsv",
        sites,
        cds_mapping,
        target=ids.target,
        uniprot=ids.uniprot,
    )

    # Extract BEB probability table from rst
    if rst.exists():
        beb_txt = workdir / f"{prefix}_beb.txt"
        try:
            extract_beb_table(rst, beb_txt)
        except RuntimeError as e:
            log.warning("BEB table extraction failed: %s", e)
    else:
        log.warning("%s missing — skipping BEB table extraction.", rst)


# =====================================================================
#  Full pipeline (standalone, no Nextflow)
# =====================================================================


def run_codeml_pipeline(
    gene_symbol: str,
    outdir: Path | None,
    skip_codeml: bool = False,
    uniprot: str | None = None,
) -> None:
    """Run both CodeML models sequentially, then extract BEB results."""
    ids = resolve_gene(gene_symbol, outdir=outdir, uniprot_override=uniprot)
    ids.workdir.mkdir(parents=True, exist_ok=True)

    prefix = ids.prefix
    workdir = ids.workdir

    phy = workdir / f"{prefix}_subset.cds.mafft.trimal.phy"
    tree = workdir / f"{prefix}_subset.tree"

    if not phy.exists() or not tree.exists():
        raise RuntimeError(
            f"Input files missing: {phy}, {tree}. Run prepare_data.py (Step 1) first."
        )

    # Write and run both models sequentially
    for model_name, model in MODELS.items():
        ctl = workdir / f"{prefix}_{model_name}.ctl"
        mlc_name = f"{prefix}_{model_name}.mlc"

        write_ctl(
            ctl,
            seqfile=phy.name,
            treefile=tree.name,
            outfile=mlc_name,
            nssites=model.nssites,
            fix_omega=model.fix_omega,
        )

        if skip_codeml:
            continue

        mlc_path = workdir / mlc_name
        run_codeml(ctl, workdir, mlc=mlc_path)

    if skip_codeml:
        log.info("--skip-codeml set – stopping after CTL generation.")
        return

    # Preserve M8 rst
    rst_src = workdir / f"{prefix}_M8_rst"
    rst_dst = workdir / f"{prefix}_M8.rst.txt"
    if rst_src.exists() and (
        not rst_dst.exists() or rst_src.stat().st_mtime >= rst_dst.stat().st_mtime
    ):
        shutil.move(str(rst_src), str(rst_dst))
        log.info("Moved %s → %s", rst_src.name, rst_dst.name)
    elif not rst_dst.exists():
        log.warning("No rst file produced by CodeML – BEB analysis will be unavailable.")

    # Extract results
    run_extract(gene_symbol, outdir=outdir, uniprot=uniprot)


# =====================================================================
#  CLI
# =====================================================================


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Step 2: Run CodeML site models (M8/M8a), report LRT, "
        "and extract BEB sites with protein positions.",
    )
    add_common_args(ap)

    mode = ap.add_mutually_exclusive_group()
    mode.add_argument(
        "--model",
        choices=list(MODELS),
        default=None,
        help="Run a single model (M8 or M8a). Designed for Nextflow "
        "parallelism. Omit to run both models sequentially.",
    )
    mode.add_argument(
        "--extract-only",
        action="store_true",
        help="Skip codeml execution; only extract BEB results from "
        "existing .mlc and .rst.txt files.",
    )

    ap.add_argument(
        "--skip-codeml",
        action="store_true",
        help="Only generate CTL files; do not run codeml. "
        "Ignored when --model or --extract-only is set.",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()

    if args.extract_only:
        run_extract(
            gene_symbol=args.gene_symbol,
            outdir=args.outdir,
            uniprot=args.uniprot,
        )
    elif args.model:
        run_single_model(
            gene_symbol=args.gene_symbol,
            model_name=args.model,
            outdir=args.outdir,
            uniprot=args.uniprot,
        )
    else:
        run_codeml_pipeline(
            gene_symbol=args.gene_symbol,
            outdir=args.outdir,
            skip_codeml=args.skip_codeml,
            uniprot=args.uniprot,
        )


if __name__ == "__main__":
    main()
