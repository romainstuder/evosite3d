#!/usr/bin/env python3
"""Step 2-ter — Run CodeML branch-site Model A along the target lineage.

The branch-site model A (Yang & Nielsen 2002; Zhang, Nielsen & Yang 2005;
``model = 2``, ``NSsites = 2``) tests for episodic positive selection on a
designated *foreground* branch. This script applies it to **every branch on
the path from the root to the target protein's leaf**, one branch at a time.

For each such branch it

    1. writes a copy of the pruned gene tree with that branch marked ``#1``,
    2. runs the alternative model MA  (``fix_omega = 0``), and
    3. runs the null model       MA0 (``fix_omega = 1``, ``omega = 1``),

then reports the likelihood-ratio test (2ΔlnL against a 50:50 mixture of a
point mass at 0 and χ²₁ — the recommended branch-site boundary correction)
and extracts the Bayes Empirical Bayes (BEB) sites for branches that test
significant.

Assumes Step 1 (``prepare_data.py``) / Step 2a have produced the PHYLIP
alignment and pruned tree in the working directory.

Outputs (per branch ``bNN``, numbered from the root towards the target):
    - ``{prefix}_bNN.tree``            — pruned tree with the branch marked #1
    - ``{prefix}_bsMA_bNN.ctl/.mlc``   — alternative model A
    - ``{prefix}_bsMA0_bNN.ctl/.mlc``  — null model A0
    - ``{prefix}_bsMA_bNN_beb_sites.tsv`` — BEB sites (significant branches)
And once per run:
    - ``{prefix}_branch_site_lrt.tsv`` — LRT summary, one row per branch

Note: the pruned Compara tree is rooted; the branch nearest the root is
labelled as given. CodeML treats the tree as unrooted (``clock = 0``), so a
basal branch may be confounded with its sister — interpret root-adjacent
branches with that caveat.

External tools required on ``PATH``:
    - ``codeml`` (PAML)

Examples::

    run_codeml_branch_site_models.py --gene-symbol HLA-DQB1
    run_codeml_branch_site_models.py --gene-symbol HLA-DQB1 --list-branches
    run_codeml_branch_site_models.py --gene-symbol HLA-DQB1 --branch 3
"""

from __future__ import annotations

import argparse
import io
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from _common import (  # type: ignore
    add_common_args,
    log,
    map_sites_to_cds,
    parse_beb_sites,
    parse_lnl,
    resolve_gene,
)
from Bio import Phylo
from run_codeml_site_models import run_codeml, write_beb_sites_tsv  # type: ignore
from scipy.stats import chi2

BRANCH_SITE_OMEGA = 1.5
PVALUE_THRESHOLD = 0.05
LRT_DF = 1


# =====================================================================
#  Branch-site control-file template
# =====================================================================


_CTL_BRANCH_SITE = """\
     seqfile = {seqfile}  * sequence data file name
    treefile = {treefile}  * tree structure file name
     outfile = {outfile}  * main result file name

       noisy = 9   * 0,1,2,3,9: how much rubbish on the screen
     verbose = 1   * 1: detailed output, 0: concise output
     runmode = 0   * 0: user tree

     seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
   CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
       clock = 0   * 0: no clock, unrooted tree
      aaDist = 0   * 0:equal, +:geometric; -:linear
       model = 2   * branch models: 2 = 2 or more dN/dS ratios for branches
     NSsites = 2   * branch-site model A (used together with model = 2)
       icode = 0   * 0:standard genetic code
       Mgene = 0   * 0:rates
   fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
       kappa = 2   * initial or fixed kappa
   fix_omega = {fix_omega}   * 1: omega fixed (null A0), 0: estimate (alt. A)
       omega = {omega}   * initial or fixed omega
       getSE = 0       * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0       * (0,1,2): rates or ancestral states
  Small_Diff = .45e-6  * Default value.
   cleandata = 0       * remove sites with ambiguity data (1:yes, 0:no)?
 fix_blength = 0       * 0: ignore branch lengths in the tree file
"""


def write_branch_site_ctl(
    ctl_path: Path,
    seqfile: str,
    treefile: str,
    outfile: str,
    fix_omega: int,
    omega: float,
) -> None:
    """Write a CodeML control file for branch-site model A (or its null)."""
    ctl_path.write_text(
        _CTL_BRANCH_SITE.format(
            seqfile=seqfile,
            treefile=treefile,
            outfile=outfile,
            fix_omega=fix_omega,
            omega=omega,
        )
    )
    log.info("Wrote %s", ctl_path)


# =====================================================================
#  Tree handling — foreground labelling along the target lineage
# =====================================================================


def load_tree(tree_path: Path) -> Any:
    """Read a Newick tree, skipping an optional PAML ``<n> 1`` header line."""
    lines = tree_path.read_text().strip().splitlines()
    if lines and re.fullmatch(r"\s*\d+\s+\d+\s*", lines[0]):
        lines = lines[1:]
    return Phylo.read(io.StringIO("".join(lines)), "newick")


def _to_newick(clade: Any, foreground: Any) -> str:
    """Serialise ``clade`` to Newick, marking ``foreground`` with ``#1``.

    Branch lengths are intentionally dropped: codeml re-estimates them
    (``fix_blength = 0``), and omitting them avoids any label/length
    ordering ambiguity in the PAML tree format.
    """
    if clade.clades:
        body = "(" + ",".join(_to_newick(c, foreground) for c in clade.clades) + ")"
    else:
        body = str(clade.name)
    if clade is foreground:
        body += " #1"
    return body


def labelled_newick(tree: Any, foreground: Any) -> str:
    """Return the full Newick string with exactly one ``#1`` foreground branch."""
    return _to_newick(tree.root, foreground) + ";"


@dataclass
class Branch:
    """One branch on the root→target path, plus its file-safe label."""

    index: int  # 1-based, counted from the root towards the target
    label: str  # filename-safe, e.g. "b03"
    description: str
    is_terminal: bool
    clade: Any  # the Bio.Phylo clade whose incoming branch is the foreground


def lineage_branches(tree: Any, target_name: str) -> list[Branch]:
    """Return every branch on the path from the root to ``target_name``."""
    target = next((t for t in tree.get_terminals() if t.name == target_name), None)
    if target is None:
        raise RuntimeError(
            f"Target leaf {target_name!r} not found in tree "
            f"({tree.count_terminals()} leaves). Check --target or the resolver."
        )
    branches: list[Branch] = []
    for i, clade in enumerate(tree.get_path(target), start=1):
        if clade.is_terminal():
            desc = f"terminal branch to {clade.name}"
        else:
            desc = f"internal branch above MRCA of {clade.count_terminals()} tips"
        branches.append(Branch(i, f"b{i:02d}", desc, clade.is_terminal(), clade))
    return branches


# =====================================================================
#  Likelihood-ratio test
# =====================================================================


def branch_site_lrt(lnl_ma: float, lnl_ma0: float) -> tuple[float, float]:
    """Return ``(2ΔlnL, p_value)`` for MA vs MA0.

    The null distribution of the branch-site test is a 50:50 mixture of a
    point mass at 0 and χ²₁, so the p-value is ``0.5 · P(χ²₁ > 2ΔlnL)``.
    """
    stat = 2.0 * (lnl_ma - lnl_ma0)
    if stat <= 0:
        return stat, 1.0
    return stat, 0.5 * float(chi2.sf(stat, LRT_DF))


@dataclass
class BranchResult:
    """Parsed LRT outcome for one branch."""

    branch: Branch
    np_ma: int
    lnl_ma: float
    np_ma0: int
    lnl_ma0: float
    stat: float
    pvalue: float
    n_beb: int
    significant: bool


# =====================================================================
#  Per-branch runner
# =====================================================================


def run_one_branch(
    branch: Branch,
    tree: Any,
    workdir: Path,
    prefix: str,
    phy: Path,
    n_leaves: int,
    omega: float,
    skip_codeml: bool,
) -> tuple[Path, Path]:
    """Write the labelled tree + both CTLs and (unless skipping) run codeml.

    Returns the ``(MA, MA0)`` ``.mlc`` paths.
    """
    tree_file = workdir / f"{prefix}_{branch.label}.tree"
    tree_file.write_text(f"{n_leaves} 1\n{labelled_newick(tree, branch.clade)}\n")

    ctl_ma = workdir / f"{prefix}_bsMA_{branch.label}.ctl"
    mlc_ma = workdir / f"{prefix}_bsMA_{branch.label}.mlc"
    write_branch_site_ctl(ctl_ma, phy.name, tree_file.name, mlc_ma.name, fix_omega=0, omega=omega)

    ctl_ma0 = workdir / f"{prefix}_bsMA0_{branch.label}.ctl"
    mlc_ma0 = workdir / f"{prefix}_bsMA0_{branch.label}.mlc"
    write_branch_site_ctl(ctl_ma0, phy.name, tree_file.name, mlc_ma0.name, fix_omega=1, omega=1)

    if skip_codeml:
        log.info("--skip-codeml: wrote CTLs + tree for %s only.", branch.label)
        return mlc_ma, mlc_ma0

    run_codeml(ctl_ma, workdir, mlc=mlc_ma)
    run_codeml(ctl_ma0, workdir, mlc=mlc_ma0)
    return mlc_ma, mlc_ma0


def _extract_branch_result(
    branch: Branch,
    mlc_ma: Path,
    mlc_ma0: Path,
    workdir: Path,
    prefix: str,
    target: str,
    uniprot: str | None,
    pvalue_threshold: float,
) -> BranchResult | None:
    """Parse the MA/MA0 mlc pair, run the LRT and (if significant) BEB sites."""
    ma = parse_lnl(mlc_ma)
    ma0 = parse_lnl(mlc_ma0)
    if not ma or not ma0:
        log.warning("Missing lnL for branch %s — skipping LRT.", branch.label)
        return None

    lnl_ma, lnl_ma0 = ma[-1], ma0[-1]
    stat, pvalue = branch_site_lrt(lnl_ma.lnl, lnl_ma0.lnl)
    significant = pvalue < pvalue_threshold

    n_beb = 0
    if significant:
        sites = parse_beb_sites(mlc_ma)
        n_beb = len(sites)
        full_aln = workdir / f"{prefix}_subset.cds.mafft.fasta"
        trimal_cols = workdir / f"{prefix}_subset.cds.mafft.trimal.cols"
        if sites and full_aln.exists() and trimal_cols.exists():
            mapping = map_sites_to_cds(full_aln, trimal_cols, target, [s.trimal_pos for s in sites])
            write_beb_sites_tsv(
                workdir / f"{prefix}_bsMA_{branch.label}_beb_sites.tsv",
                sites,
                mapping,
                target=target,
                uniprot=uniprot,
            )
        elif sites:
            log.warning(
                "Cannot map trimal→CDS positions for %s (missing %s or %s).",
                branch.label,
                full_aln,
                trimal_cols,
            )

    return BranchResult(
        branch=branch,
        np_ma=lnl_ma.np,
        lnl_ma=lnl_ma.lnl,
        np_ma0=lnl_ma0.np,
        lnl_ma0=lnl_ma0.lnl,
        stat=stat,
        pvalue=pvalue,
        n_beb=n_beb,
        significant=significant,
    )


# =====================================================================
#  Reporting
# =====================================================================


def write_lrt_summary(path: Path, results: list[BranchResult], target: str) -> None:
    """Write the per-branch LRT summary TSV."""
    lines = [
        f"# Branch-site model A (MA vs MA0) LRT along the lineage to {target}.",
        "# p_value = 0.5 * P(chi2_1 > 2DeltaLnL) (50:50 boundary mixture). df = 1.",
        "branch\tdescription\tnp_MA\tlnL_MA\tnp_MA0\tlnL_MA0\t"
        "2DeltaLnL\tdf\tp_value\tn_BEB_sites\tsignificant",
    ]
    for r in results:
        lines.append(
            f"{r.branch.label}\t{r.branch.description}\t{r.np_ma}\t{r.lnl_ma:.4f}\t"
            f"{r.np_ma0}\t{r.lnl_ma0:.4f}\t{r.stat:.4f}\t{LRT_DF}\t{r.pvalue:.3e}\t"
            f"{r.n_beb}\t{'yes' if r.significant else 'no'}"
        )
    path.write_text("\n".join(lines) + "\n")
    log.info("Wrote LRT summary (%d branches) → %s", len(results), path)


def print_report(results: list[BranchResult], target: str) -> None:
    """Print a console summary of the branch-site LRTs."""
    print(f"\n=== Branch-site Model A LRT — lineage to {target} ===")
    print(f"{'branch':>6} {'2dLnL':>10} {'p-value':>10} {'BEB':>4}  description")
    for r in results:
        flag = "*" if r.significant else " "
        print(
            f"{r.branch.label:>6} {r.stat:10.4f} {r.pvalue:10.3e} "
            f"{r.n_beb:>4}{flag} {r.branch.description}"
        )
    n_sig = sum(1 for r in results if r.significant)
    print(f"\n{n_sig} / {len(results)} branches significant.\n")


# =====================================================================
#  Pipeline
# =====================================================================


def run_branch_site_pipeline(
    gene_symbol: str,
    outdir: Path | None,
    uniprot: str | None = None,
    target: str | None = None,
    tree_path: Path | None = None,
    omega: float = BRANCH_SITE_OMEGA,
    only_branch: int | None = None,
    skip_codeml: bool = False,
    pvalue_threshold: float = PVALUE_THRESHOLD,
) -> None:
    """Run branch-site model A on every branch from the root to the target."""
    ids = resolve_gene(gene_symbol, outdir=outdir, uniprot_override=uniprot)
    prefix, workdir = ids.prefix, ids.workdir
    workdir.mkdir(parents=True, exist_ok=True)

    phy = workdir / f"{prefix}_subset.cds.mafft.trimal.phy"
    tree_file_in = tree_path if tree_path is not None else workdir / f"{prefix}_subset.tree"
    if not phy.exists() or not tree_file_in.exists():
        raise RuntimeError(
            f"Input files missing: {phy}, {tree_file_in}. "
            "Run prepare_data.py (Step 1) and the Silver transform (Step 2a) first."
        )

    target_name = target or ids.target
    tree = load_tree(tree_file_in)
    n_leaves = tree.count_terminals()
    branches = lineage_branches(tree, target_name)
    log.info("Lineage to %s has %d branches from the root.", target_name, len(branches))

    if only_branch is not None:
        branches = [b for b in branches if b.index == only_branch]
        if not branches:
            raise RuntimeError(f"--branch {only_branch} is out of range (1..{n_leaves}).")

    results: list[BranchResult] = []
    for branch in branches:
        log.info("=== Branch %s: %s ===", branch.label, branch.description)
        mlc_ma, mlc_ma0 = run_one_branch(
            branch, tree, workdir, prefix, phy, n_leaves, omega, skip_codeml
        )
        if skip_codeml:
            continue
        result = _extract_branch_result(
            branch, mlc_ma, mlc_ma0, workdir, prefix, target_name, ids.uniprot, pvalue_threshold
        )
        if result is not None:
            results.append(result)

    if skip_codeml:
        log.info("--skip-codeml set – stopping after CTL/tree generation.")
        return

    if results:
        write_lrt_summary(workdir / f"{prefix}_branch_site_lrt.tsv", results, target_name)
        print_report(results, target_name)


def list_branches(
    gene_symbol: str,
    outdir: Path | None,
    uniprot: str | None,
    target: str | None,
    tree_path: Path | None,
) -> None:
    """Print the lineage branches for a gene and return without running codeml."""
    ids = resolve_gene(gene_symbol, outdir=outdir, uniprot_override=uniprot)
    tree_file_in = tree_path if tree_path is not None else ids.workdir / f"{ids.prefix}_subset.tree"
    if not tree_file_in.exists():
        raise RuntimeError(
            f"Tree not found: {tree_file_in}. Run the Silver transform (Step 2a) first."
        )
    tree = load_tree(tree_file_in)
    target_name = target or ids.target
    branches = lineage_branches(tree, target_name)
    print(f"\nLineage from root to {target_name} — {len(branches)} branches:")
    for b in branches:
        print(f"  {b.index:>3}  {b.label}  {b.description}")
    print()


# =====================================================================
#  CLI
# =====================================================================


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Step 2-ter: Run CodeML branch-site model A (MA vs MA0) on "
        "every branch of the phylogeny from the root to the target protein.",
    )
    add_common_args(ap)
    ap.add_argument(
        "--target",
        default=None,
        help="Ensembl protein ID of the foreground leaf (default: the resolver's protein).",
    )
    ap.add_argument(
        "--tree",
        type=Path,
        default=None,
        help="Override the pruned tree (default: {prefix}_subset.tree in the workdir).",
    )
    ap.add_argument(
        "--omega",
        type=float,
        default=BRANCH_SITE_OMEGA,
        help=f"Initial omega for the alternative model (default {BRANCH_SITE_OMEGA}).",
    )
    ap.add_argument(
        "--branch",
        type=int,
        default=None,
        help="Run only this 1-based branch index (root→target). Default: all branches.",
    )
    ap.add_argument(
        "--pvalue-threshold",
        type=float,
        default=PVALUE_THRESHOLD,
        help=f"LRT significance threshold for flagging / BEB extraction (default {PVALUE_THRESHOLD}).",
    )
    ap.add_argument(
        "--list-branches",
        action="store_true",
        help="List the lineage branches and exit (no codeml).",
    )
    ap.add_argument(
        "--skip-codeml",
        action="store_true",
        help="Write CTL files and labelled trees but do not run codeml.",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    if args.list_branches:
        list_branches(
            gene_symbol=args.gene_symbol,
            outdir=args.outdir,
            uniprot=args.uniprot,
            target=args.target,
            tree_path=args.tree,
        )
        return
    run_branch_site_pipeline(
        gene_symbol=args.gene_symbol,
        outdir=args.outdir,
        uniprot=args.uniprot,
        target=args.target,
        tree_path=args.tree,
        omega=args.omega,
        only_branch=args.branch,
        skip_codeml=args.skip_codeml,
        pvalue_threshold=args.pvalue_threshold,
    )


if __name__ == "__main__":
    main()
