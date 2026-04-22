#!/usr/bin/env python3
"""Run CodeML site-models for pervasive positive selection.

Takes a plain HGNC gene symbol (e.g. ``HLA-DQB1``) and resolves it to
the Ensembl gene, canonical protein stable ID and UniProt accession via
the shared resolver in
``software/protein_characteriser/bin/resolve_ids.py``, then automates
part 2 of the tutorial ``tutorials/pervasive_selection_site_models``:

    1. Ensuring the CDS PHYLIP alignment and pruned gene tree are
       present (calls :mod:`scripts.fetch_ensembl_msa_tree` if not).
       The fetch script filters to the curated reference-species panel
       (41 vertebrates, fishes to mammals) *before* TrimAl, so gap
       statistics reflect only the species used for CodeML.
    2. Writing two CodeML control files:
           - ``{prefix}_M8.ctl``   (NSsites = 8)
           - ``{prefix}_M8a.ctl``  (NSsites = 8, fix_omega = 1)
    4. Running ``codeml`` on each control file.
    5. Preserving ``rst`` as ``{prefix}_M8.rst.txt``.
    6. Reporting lnL per model and the M8a-M8 likelihood-ratio test
       (the preferred site-model test for pervasive positive selection).
    7. Post-analysis (skip with ``--skip-post-analysis``):
           - Extract BEB ``Positively selected sites`` from the M8
             ``.mlc`` and map each trimmed column back to its CDS /
             reference-protein position via
             :mod:`scripts.get_position_cds_trimal`
             (``{prefix}_beb_sites.tsv``).
           - Extract the per-site BEB probability block from the rst
             file, normalise its spacing and run
             :mod:`scripts.plot_dnds_per_site`
             (``{prefix}_beb.txt``, ``{prefix}_beb.png``).
           - Write a Jalview annotation file with BEB Pr(w>1) bar
             graph and posterior mean dN/dS line graph tracks
             (``{prefix}_beb.jlv``).
           - Fetch an AlphaFold monomer (or an RCSB PDB via ``--pdb``)
             for the target and write ``{prefix}_sites.pml`` — a PyMOL
             script that colours BEB-significant residues as yellow
             spheres / sticks.

External tools required on ``PATH``:
    - ``codeml`` (PAML)
    - ``nw_prune`` (Newick Utilities) — used to re-prune the gene tree
      after the reference-species filter.
    - ``trimal`` (only needed if the inputs are missing and
      :mod:`scripts.fetch_ensembl_msa_tree` has to run).

Network access to ``rest.ensembl.org``, ``alphafold.ebi.ac.uk`` and
(optionally) ``files.rcsb.org`` is required for the post-analysis step.

Example:
    run_codeml_site_models.py HLA-DQB1
    run_codeml_site_models.py HLA-DQB1 --outdir HLA_DQB1_workdir --uniprot P01920
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import logging
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from scipy.stats import chi2

_REPO = Path(__file__).resolve().parent.parent.parent
_FETCH_SCRIPT = _REPO / "scripts" / "fetch_ensembl_msa_tree.py"
_GET_POS_SCRIPT = _REPO / "scripts" / "get_position_cds_trimal.py"
_PLOT_SCRIPT = _REPO / "scripts" / "plot_dnds_per_site.py"

# Reuse the gene-symbol → Ensembl / UniProt resolver from
# software/protein_characteriser so a single entry point can take a
# plain gene symbol (e.g. ``HLA-DQB1``) and derive the Ensembl gene,
# canonical protein stable ID and UniProt accession.
_PC = _REPO
sys.path.insert(0, str(_PC / "bin"))
sys.path.insert(0, str(_PC / "lib"))
sys.path.insert(0, str(_PC / "scripts"))

from resolve_ids import resolve as resolve_identifier  # noqa: E402

_ALPHAFOLD_FILES = "https://alphafold.ebi.ac.uk/files"
_ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"
_RCSB_FILES = "https://files.rcsb.org/download"

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


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
    """Write a CodeML control file.

    Args:
        ctl_path: Destination for the control file.
        seqfile: Path (relative to the CTL) to the PHYLIP alignment.
        treefile: Path (relative to the CTL) to the Newick tree.
        outfile: Name of the ``.mlc`` output file.
        nssites: Space-separated list of site-model indices (e.g. ``"8"``).
        fix_omega: ``1`` to fix omega (M8a null), ``0`` otherwise.
    """
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
#  Input preparation
# =====================================================================


def _check_tool(name: str) -> None:
    """Raise ``RuntimeError`` if an external binary is missing from ``PATH``."""
    if shutil.which(name) is None:
        raise RuntimeError(f"Required tool '{name}' not found on PATH.")


def ensure_inputs(
    gene: str,
    prefix: str,
    target: str | None,
    outdir: Path,
    taxon: int,
    phy: Path,
    tree: Path,
) -> None:
    """Make sure the PHYLIP alignment and pruned tree exist.

    Invokes :mod:`scripts.fetch_ensembl_msa_tree` when either file is
    missing. No-op when both already exist.

    Args:
        gene: Ensembl gene stable ID.
        prefix: Filename stem shared with the fetch script.
        target: Ensembl protein ID to place first in the MSA.
        outdir: Working directory for all outputs.
        taxon: NCBI taxon ID to prune the Compara gene tree.
        phy: Expected PHYLIP alignment path.
        tree: Expected pruned-tree path.
    """
    if phy.exists() and tree.exists():
        log.info("Found existing inputs – skipping fetch.")
        return

    log.info("Missing inputs – running fetch_ensembl_msa_tree.py")
    cmd = [
        sys.executable,
        str(_FETCH_SCRIPT),
        "--gene",
        gene,
        "--prefix",
        prefix,
        "--outdir",
        str(outdir),
        "--taxon",
        str(taxon),
    ]
    if target:
        cmd += ["--target", target]
    subprocess.run(cmd, check=True)

    if not phy.exists() or not tree.exists():
        raise RuntimeError(
            f"fetch_ensembl_msa_tree.py did not produce expected files: {phy}, {tree}"
        )


# =====================================================================
#  Running CodeML
# =====================================================================


def run_codeml(ctl: Path, workdir: Path, mlc: Path | None = None) -> None:
    """Run ``codeml`` on a control file in an isolated subdirectory.

    CodeML writes fixed-name scratch files (``rst``, ``rst1``, ``rub``,
    ``2NG.dN``, …) into the working directory, so concurrent runs must
    each get their own folder. This function creates a subdirectory
    named after the CTL stem (e.g. ``HLA_DQB1_M8/``), symlinks the
    input alignment and tree, copies the CTL file, runs ``codeml``
    there, and moves the outputs back into ``workdir``.

    Args:
        ctl: Path to the control file (must live inside ``workdir``).
        workdir: Parent directory that holds the shared inputs.
        mlc: Expected ``.mlc`` output in ``workdir``. When provided and
            already present, the CodeML invocation is skipped.
    """
    if mlc is not None and mlc.exists():
        log.info("Skipping codeml on %s – %s already exists.", ctl.name, mlc.name)
        return
    _check_tool("codeml")

    # Create an isolated run directory so scratch files don't collide.
    rundir = workdir / ctl.stem
    rundir.mkdir(exist_ok=True)

    # Symlink shared inputs (alignment + tree) referenced by the CTL.
    for src in workdir.iterdir():
        if src.suffix in (".phy", ".tree") and src.is_file():
            dst = rundir / src.name
            if not dst.exists():
                dst.symlink_to(src.resolve())

    # Copy the CTL file into the run directory.
    shutil.copy2(ctl, rundir / ctl.name)

    log.info("Running codeml on %s (this can take some minutes)", ctl.name)
    subprocess.run(["codeml", ctl.name], check=True, cwd=rundir)

    # Move outputs back into the main workdir.  Scratch files with
    # fixed names (rst, rst1, rub, 2NG.dN, …) are prefixed with the
    # CTL stem so parallel runs don't overwrite each other.
    mlc_names = {f for f in rundir.iterdir() if f.suffix == ".mlc"}
    for path in sorted(rundir.iterdir()):
        if path.is_symlink() or path.name == ctl.name:
            continue
        if path in mlc_names:
            dest = workdir / path.name
        else:
            dest = workdir / f"{ctl.stem}_{path.name}"
        shutil.move(str(path), str(dest))

    # Clean up symlinks left in the run directory.
    for path in rundir.iterdir():
        path.unlink()
    rundir.rmdir()


# =====================================================================
#  Parsing lnL and LRTs
# =====================================================================


_LNL_RE = re.compile(r"lnL\(ntime:\s*(\d+)\s+np:\s*(\d+)\):\s*(-?\d+\.\d+)")


@dataclass
class LnL:
    """A single ``lnL`` line extracted from an ``.mlc`` file."""

    ntime: int
    np: int
    lnl: float


def parse_lnl(mlc: Path) -> list[LnL]:
    """Extract every ``lnL`` line from a CodeML ``.mlc`` file, in order.

    Args:
        mlc: Path to the CodeML output file.

    Returns:
        List of :class:`LnL` entries, one per site-model computed.
    """
    entries = []
    for line in mlc.read_text().splitlines():
        m = _LNL_RE.search(line)
        if m:
            entries.append(LnL(int(m.group(1)), int(m.group(2)), float(m.group(3))))
    return entries


def lrt(name: str, null: LnL, alt: LnL) -> str:
    """Format a single likelihood-ratio test line.

    Args:
        name: Human-readable LRT name (e.g. ``"M1a-M2a"``).
        null: Null-model lnL entry.
        alt: Alternative-model lnL entry.

    Returns:
        Formatted summary string with 2ΔlnL, df and p-value.
    """
    stat = 2.0 * (alt.lnl - null.lnl)
    df = alt.np - null.np
    p = chi2.sf(stat, df) if df > 0 and stat > 0 else float("nan")
    return f"{name:>7}  2ΔlnL = {stat:10.4f}   df = {df}   p = {p:.3e}"


def report_lrts(prefix: str, outdir: Path) -> None:
    """Parse the CodeML outputs and print the M8a-M8 LRT.

    Args:
        prefix: Filename stem used for CodeML outputs.
        outdir: Directory containing the ``.mlc`` files.
    """
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
#  BEB site extraction
# =====================================================================


@dataclass
class BebSite:
    """A positively-selected site from a CodeML BEB table.

    Attributes:
        trimal_pos: 1-based position in the trimmed (CodeML) alignment.
        aa: Reference-sequence amino acid at that column.
        prob: ``Pr(w > 1)`` posterior probability.
        stars: Significance marker (``""``, ``"*"`` for P > 95 %, ``"**"``
            for P > 99 %).
        mean_w: Posterior-mean dN/dS for that site.
        se: Posterior standard error of ``mean_w``.
    """

    trimal_pos: int
    aa: str
    prob: float
    stars: str
    mean_w: float
    se: float


_BEB_HEADER = "Bayes Empirical Bayes (BEB) analysis"
_BEB_SITE_RE = re.compile(
    r"^\s*(\d+)\s+([A-Z*\-])\s+(\d+\.\d+)(\*{0,2})\s+"
    r"(-?\d+\.\d+)\s+\+-\s+(-?\d+\.\d+)\s*$"
)


def parse_beb_sites(mlc: Path) -> list[BebSite]:
    """Extract the BEB ``Positively selected sites`` block from an ``.mlc`` file.

    Scans for the *last* ``Bayes Empirical Bayes (BEB) analysis`` header
    — an ``M8`` run has a single such block; a combined M2a/M8 run has
    two and the M8 one is always last. Reads subsequent ``<pos> <aa>
    <prob>[*|**] <mean> +- <se>`` rows until the block terminates on a
    blank line.

    Args:
        mlc: Path to the CodeML output file.

    Returns:
        List of :class:`BebSite` entries in file order. Empty if no BEB
        block is found.
    """
    text = mlc.read_text()
    idx = text.rfind(_BEB_HEADER)
    if idx < 0:
        log.warning("No BEB block found in %s", mlc)
        return []

    sites: list[BebSite] = []
    in_table = False
    for line in text[idx:].splitlines()[1:]:
        m = _BEB_SITE_RE.match(line)
        if m:
            in_table = True
            sites.append(
                BebSite(
                    trimal_pos=int(m.group(1)),
                    aa=m.group(2),
                    prob=float(m.group(3)),
                    stars=m.group(4),
                    mean_w=float(m.group(5)),
                    se=float(m.group(6)),
                )
            )
        elif in_table and not line.strip():
            break
    return sites


# =====================================================================
#  Mapping BEB sites back to CDS / reference-protein positions
# =====================================================================


_GET_POS_RE = re.compile(r"site_index:\s*(\d+)\s+CDS pos:\s*(\d+|N/A)\s+AA:\s*(\S+)")


def map_sites_to_cds(
    full_aln: Path,
    trimal_cols: Path,
    target: str,
    sites: list[int],
) -> dict[int, tuple[int | None, str]]:
    """Map CodeML/trimal column indices to CDS positions in the target.

    Delegates to :mod:`scripts.get_position_cds_trimal` so the tutorial's
    column-indexing conventions are preserved.

    Args:
        full_aln: Pre-TrimAl CDS alignment (``{prefix}_subset.cds.mafft.fasta``).
        trimal_cols: TrimAl ``-colnumbering`` output
            (``{prefix}_subset.cds.mafft.trimal.cols``).
        target: Reference sequence ID (e.g. ``"ENSP00000407332"``).
        sites: Trimmed-alignment site indices from the BEB table.

    Returns:
        Mapping ``{trimal_pos: (cds_pos, aa)}``. ``cds_pos`` is ``None``
        when the helper script could not resolve a position.
    """
    if not sites:
        return {}
    site_arg = " ".join(str(s) for s in sites)
    cmd = [
        sys.executable,
        str(_GET_POS_SCRIPT),
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
    """Write a tab-separated summary of BEB sites + target-protein positions.

    Args:
        path: Destination TSV path.
        sites: BEB sites parsed from the ``.mlc`` file.
        cds_mapping: Output of :func:`map_sites_to_cds`.
        target: Ensembl protein stable ID of the reference sequence (used
            only as metadata in the TSV header).
        uniprot: UniProt accession of the reference sequence. When the
            Ensembl translation and the UniProt canonical isoform match
            residue-for-residue (the usual case for canonical protein
            entries), ``protein_pos`` equals the UniProt position too.

    The ``protein_pos`` column is 1-based along the ungapped target
    protein sequence — i.e. the residue numbering shared by the Ensembl
    translation and the UniProt canonical sequence — making the file
    directly usable against either ID.
    """
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
#  Per-site dN/dS plot (Manhattan-style)
# =====================================================================


_BEB_TABLE_HEADER_RE = re.compile(
    r"(Bayes Empirical Bayes \(BEB\)|Naive Empirical Bayes \(NEB\))"
    r" probabilities for (\d+) classes"
)
_BEB_TABLE_ROW_RE = re.compile(r"^\s*\d+\s+\S\s+")


def extract_beb_table(rst_txt: Path, out_beb_txt: Path) -> int:
    """Extract the per-site posterior-probabilities block from ``rst``.

    Locates the ``Bayes Empirical Bayes (BEB) probabilities for N
    classes`` table in the rst file, keeps only the data rows, and
    normalises the ``( N)`` → ``(N)`` spacing so every row has the same
    whitespace-separated column count (required by
    :mod:`scripts.plot_dnds_per_site`).

    BEB is preferred, but a ``Naive Empirical Bayes (NEB)`` block is
    used as fallback — some PAML runs (notably those with ``ncatG``
    reduced from the default ``10`` to ``5``) produce a per-site NEB
    table but no per-site BEB table in the rst.

    ``plot_dnds_per_site.py`` hardcodes ``iloc[:, 12]`` for the
    positive-selection class probability and ``iloc[:, 14]`` for the
    posterior mean ω, which matches the 11-class BEB layout. When the
    block found has fewer than 11 classes, each row is left-padded with
    ``0.0`` columns so the positive-selection class still lands at
    col 12.

    Args:
        rst_txt: Path to the preserved ``rst`` file
            (``{prefix}_M8.rst.txt``).
        out_beb_txt: Destination file for the cleaned table.

    Returns:
        The number of per-site rows written.

    Raises:
        RuntimeError: if no BEB/NEB probabilities block can be located.
    """
    text = rst_txt.read_text()

    # Walk all header matches; prefer the *last* BEB block, and fall
    # back to the last NEB block if no BEB block exists.
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
            # Collapse "( 3)" → "(3)" so single- and double-digit class
            # labels tokenise to one column.
            rows.append(line.replace("( ", "("))
        elif seen_data and not line.strip():
            break

    if not rows:
        raise RuntimeError(f"{which} probability table was empty in {rst_txt}")

    # Pad shorter tables so the positive-selection class sits at col 12,
    # matching the layout expected by plot_dnds_per_site.py.
    if n_classes < 11:
        pad_cols = 11 - n_classes
        pad = "0.0 " * pad_cols
        padded: list[str] = []
        for row in rows:
            # Split into (pos, aa, rest) so we can insert the pad
            # between the amino-acid column and the first probability.
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


def plot_beb_per_site(beb_txt: Path, out_png: Path, threshold: float = 0.50) -> None:
    """Run :mod:`scripts.plot_dnds_per_site` to produce the dN/dS plot.

    Args:
        beb_txt: Cleaned BEB table from :func:`extract_beb_table`.
        out_png: Destination image path.
        threshold: BEB probability threshold above which sites are coloured
            red in the plot (default 0.50, matching the tutorial).
    """
    cmd = [
        sys.executable,
        str(_PLOT_SCRIPT),
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
#  Jalview annotation tracks
# =====================================================================


def write_jalview_annotation(
    path: Path,
    beb_txt: Path,
    threshold: float = 0.50,
) -> None:
    """Write a Jalview annotation file with BEB and dN/dS tracks.

    Produces two annotation tracks loadable as alignment annotations in
    Jalview (``File → Load Features / Annotations``):

        - **BEB Pr(w>1)** — bar graph of the posterior probability that
          each site belongs to the positive-selection class (column 12
          in the padded BEB table). Bars above ``threshold`` are coloured
          red, the rest black.
        - **dN/dS (posterior mean)** — line graph of the posterior mean
          omega per site (column 14). A reference line at ω = 1 is added
          via a ``GRAPHLINE`` directive.

    Positions are 1-based alignment columns matching the trimmed
    (CodeML) alignment loaded as the corresponding ``.phy`` file.

    Args:
        path: Destination ``.jlv`` annotation file.
        beb_txt: Cleaned BEB table from :func:`extract_beb_table`
            (whitespace-separated, one row per site).
        threshold: BEB probability cutoff; sites above this are
            coloured red in the bar graph.
    """
    rows = beb_txt.read_text().splitlines()
    if not rows:
        log.warning("Empty BEB table — skipping Jalview annotation.")
        return

    beb_probs: list[float] = []
    mean_omegas: list[float] = []
    for row in rows:
        cols = row.split()
        # Column 12 = Pr(positive selection), column 14 = posterior mean ω.
        beb_probs.append(float(cols[12]))
        mean_omegas.append(float(cols[14]))

    # Build bar-graph values for BEB Pr(w>1).
    beb_vals: list[str] = []
    for p in beb_probs:
        colour = "ff0000" if p > threshold else "000000"
        beb_vals.append(f"{p:.4f},{p:.4f},{colour}")

    # Build line-graph values for dN/dS.
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
#  3D structure + PyMOL script
# =====================================================================


def _http_get(url: str, timeout: int = 60) -> bytes | None:
    """Simple HTTP GET returning raw bytes, or ``None`` on failure."""
    try:
        req = Request(url, headers={"User-Agent": "evosite3d/run_codeml_site_models"})
        with urlopen(req, timeout=timeout) as resp:
            return resp.read()
    except (HTTPError, URLError, TimeoutError) as e:
        log.warning("HTTP GET failed for %s: %s", url, e)
        return None


def fetch_alphafold_pdb(uniprot: str, out_path: Path) -> bool:
    """Download the AlphaFold monomer PDB for a UniProt accession.

    First queries the AlphaFold prediction API to discover the *latest*
    model URL (the ``pdbUrl`` field), since AlphaFold DB has cycled
    through versions v1 → v6 and older model files are eventually
    removed. Falls back to a v6 → v1 version sweep if the API is
    unreachable.

    Args:
        uniprot: UniProt accession.
        out_path: Destination PDB path.

    Returns:
        ``True`` on success.
    """
    pdb_url: str | None = None
    api_raw = _http_get(f"{_ALPHAFOLD_API}/{uniprot}")
    if api_raw is not None:
        try:
            data = json.loads(api_raw)
        except json.JSONDecodeError:
            data = None
        if isinstance(data, list) and data:
            pdb_url = data[0].get("pdbUrl")

    if pdb_url:
        raw = _http_get(pdb_url)
        if raw is not None:
            out_path.write_bytes(raw)
            log.info("AlphaFold model %s → %s", uniprot, out_path)
            return True

    # Fallback sweep from newest to oldest known version.
    for version in range(6, 0, -1):
        url = f"{_ALPHAFOLD_FILES}/AF-{uniprot}-F1-model_v{version}.pdb"
        raw = _http_get(url)
        if raw is not None:
            out_path.write_bytes(raw)
            log.info("AlphaFold model %s (v%d) → %s", uniprot, version, out_path)
            return True
    return False


def fetch_rcsb_pdb(pdb_id: str, out_path: Path) -> bool:
    """Download an experimental PDB structure from the RCSB.

    Args:
        pdb_id: 4-character PDB ID (case-insensitive).
        out_path: Destination PDB path.

    Returns:
        ``True`` on success.
    """
    url = f"{_RCSB_FILES}/{pdb_id.upper()}.pdb"
    raw = _http_get(url)
    if raw is None:
        return False
    out_path.write_bytes(raw)
    log.info("RCSB PDB %s → %s", pdb_id.upper(), out_path)
    return True


def write_pymol_script(
    pml: Path,
    pdb_path: Path,
    sites: list[BebSite],
    cds_mapping: dict[int, tuple[int | None, str]],
    chain: str = "A",
    resi_offset: int = 0,
    threshold_high: float = 0.95,
    threshold_low: float = 0.50,
) -> None:
    """Write a PyMOL ``.pml`` that highlights BEB-positive sites on the structure.

    Follows the tutorial's colour scheme:

        - yellow spheres for sites with ``Pr(w > 1) >= threshold_high``
        - yellow sticks for ``threshold_low <= Pr(w > 1) < threshold_high``

    Residues are selected by the target-protein position (the
    ``protein_pos`` column in ``{prefix}_beb_sites.tsv``), optionally
    shifted by ``resi_offset`` to account for numbering mismatches
    between the reference protein and a crystallographic construct
    (e.g. the tutorial's ``+32`` on 1UVQ). AlphaFold models share the
    UniProt canonical numbering, so they should be loaded with
    ``resi_offset = 0``.

    Args:
        pml: Destination PyMOL script.
        pdb_path: Structure file to load (same directory as ``pml`` is
            typical so the script is portable).
        sites: BEB sites parsed from the M8 ``.mlc``.
        cds_mapping: Output of :func:`map_sites_to_cds`. Values carry
            the 1-based position along the target protein sequence and
            the reference amino acid.
        chain: Chain identifier of the target protein in the structure.
        resi_offset: Offset added to every residue number in the loaded
            structure via ``alter target, resi = str(int(resi)+N)`` so
            that crystal-style numbering (e.g. 1UVQ starting at residue 1
            of the mature β1 domain) lines up with the target-protein
            numbering (which includes the signal peptide). After the
            ``alter``/``rebuild``, the BEB selections pick residues by
            their unaltered target-protein position.
        threshold_high: Cutoff for "sphere" sites.
        threshold_low: Cutoff for "stick" sites.
    """
    high, low = [], []
    for s in sites:
        protein_pos, _ = cds_mapping.get(s.trimal_pos, (None, s.aa))
        if protein_pos is None:
            continue
        if s.prob >= threshold_high:
            high.append(protein_pos)
        elif s.prob >= threshold_low:
            low.append(protein_pos)

    high_expr = "+".join(str(r) for r in high) if high else ""
    low_expr = "+".join(str(r) for r in low) if low else ""

    # Use the absolute path so the .pml works regardless of the
    # directory PyMOL is launched from.
    pdb_abs = pdb_path.resolve()
    lines = [
        "# PyMOL visualisation of BEB-positive sites (M8 site model).",
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
    if high_expr:
        lines += [
            f"select sites_BEB{int(threshold_high * 100)}, target and resi {high_expr}",
            f"show spheres, sites_BEB{int(threshold_high * 100)}",
            f"colour yellow, sites_BEB{int(threshold_high * 100)}",
        ]
    if low_expr:
        lines += [
            f"select sites_BEB{int(threshold_low * 100)}, target and resi {low_expr}",
            f"show sticks, sites_BEB{int(threshold_low * 100)}",
            f"colour orange, sites_BEB{int(threshold_low * 100)}",
        ]
    lines += [
        "bg_color white",
        "util.cnc",
        "zoom target",
        # "set ray_opaque_background, 1",
    ]
    pml.write_text("\n".join(lines) + "\n")
    log.info(
        "Wrote PyMOL script → %s (%d sphere sites, %d stick sites)",
        pml,
        len(high),
        len(low),
    )


def fetch_structure(
    target: str,
    outdir: Path,
    prefix: str,
    pdb_id: str | None,
    uniprot: str | None = None,
) -> tuple[Path | None, str]:
    """Fetch a 3D structure file for the target, preferring AlphaFold.

    Args:
        target: Ensembl protein ID for the reference sequence.
        outdir: Destination directory.
        prefix: Filename stem for AlphaFold-derived downloads.
        pdb_id: If provided, override the AlphaFold lookup and pull this
            PDB entry from the RCSB instead.
        uniprot: UniProt accession resolved upstream by
            :func:`resolve_identifier`. Required for the AlphaFold fetch.

    Returns:
        ``(path, source)`` — the on-disk structure path (or ``None`` on
        failure) and a short source label (``"alphafold"`` or ``"rcsb"``).
    """
    if pdb_id:
        dest = outdir / f"{pdb_id.upper()}.pdb"
        got = fetch_rcsb_pdb(pdb_id, dest)
        return (dest if got else None), "rcsb"

    if not uniprot:
        log.warning(
            "Could not resolve a UniProt accession for %s — "
            "skipping AlphaFold download. Use --pdb to supply a "
            "structure manually.",
            target,
        )
        return None, "alphafold"

    dest = outdir / f"{prefix}_{uniprot}_AF.pdb"
    got = fetch_alphafold_pdb(uniprot, dest)
    return (dest if got else None), "alphafold"


# =====================================================================
#  Post-analysis orchestrator
# =====================================================================


def run_post_analysis(
    prefix: str,
    target: str | None,
    outdir: Path,
    uniprot: str | None = None,
    pdb_id: str | None = None,
    chain: str = "A",
    resi_offset: int = 0,
    plot_threshold: float = 0.50,
) -> None:
    """Run the tutorial's post-CodeML analysis steps.

    Produces:
        - ``{prefix}_beb_sites.tsv`` — BEB sites with CDS positions.
        - ``{prefix}_beb.txt``       — cleaned BEB probabilities table.
        - ``{prefix}_beb.png``       — per-site dN/dS plot.
        - ``{prefix}_beb.jlv``       — Jalview annotation (BEB + dN/dS tracks).
        - ``{prefix}_sites.pml``     — PyMOL script colouring BEB sites.
        - ``{prefix}_<uniprot>_AF.pdb`` or ``{PDB_ID}.pdb`` — structure.

    Args:
        prefix: Filename stem used throughout the pipeline.
        target: Ensembl protein ID of the reference sequence (used for
            :func:`get_position_cds_trimal` and AlphaFold lookup).
        outdir: Working directory.
        uniprot: Explicit UniProt accession; when ``None`` the accession
            is resolved via :func:`fetch_uniprot_from_ensembl`.
        pdb_id: Optional PDB accession; overrides the AlphaFold fetch.
        chain: Chain to highlight in the PyMOL script.
        resi_offset: Residue-numbering offset added inside the PyMOL script
            (e.g. ``32`` for the tutorial's 1UVQ example).
        plot_threshold: BEB probability threshold for the Manhattan plot
            (default 0.50).
    """
    if not target:
        log.warning(
            "No --target supplied — skipping post-analysis (BEB mapping "
            "and structure visualisation both need a reference sequence)."
        )
        return

    mlc = outdir / f"{prefix}_M8.mlc"
    rst = outdir / f"{prefix}_M8.rst.txt"
    full_aln = outdir / f"{prefix}_subset.cds.mafft.fasta"
    trimal_cols = outdir / f"{prefix}_subset.cds.mafft.trimal.cols"

    if not mlc.exists():
        log.warning("%s missing — cannot run post-analysis.", mlc)
        return

    # 1) BEB sites → CDS positions
    sites = parse_beb_sites(mlc)
    if not sites:
        log.warning("No BEB sites parsed from %s — skipping post-analysis.", mlc)
        return
    log.info("Parsed %d BEB sites from %s", len(sites), mlc.name)

    cds_mapping: dict[int, tuple[int | None, str]] = {}
    if full_aln.exists() and trimal_cols.exists():
        cds_mapping = map_sites_to_cds(
            full_aln,
            trimal_cols,
            target,
            [s.trimal_pos for s in sites],
        )
    else:
        log.warning(
            "Cannot map trimal→CDS positions (missing %s or %s).",
            full_aln,
            trimal_cols,
        )

    # The UniProt accession is resolved upstream in ``run_pipeline`` via
    # ``resolve_identifier``; the caller may also override it with
    # ``--uniprot`` (needed when the Ensembl cross-reference lands on an
    # isoform-level TrEMBL entry such as Q5Y7D6 for HLA-DQB1, where the
    # canonical Swiss-Prot entry is P01920).
    if uniprot:
        log.info("Using UniProt accession %s for %s", uniprot, target)

    write_beb_sites_tsv(
        outdir / f"{prefix}_beb_sites.tsv",
        sites,
        cds_mapping,
        target=target,
        uniprot=uniprot,
    )

    # 2) dN/dS plot + Jalview annotation
    if rst.exists():
        beb_txt = outdir / f"{prefix}_beb.txt"
        beb_png = outdir / f"{prefix}_beb.png"
        jlv = outdir / f"{prefix}_beb.jlv"
        try:
            extract_beb_table(rst, beb_txt)
            plot_beb_per_site(beb_txt, beb_png, threshold=plot_threshold)
            write_jalview_annotation(jlv, beb_txt, threshold=plot_threshold)
        except (RuntimeError, subprocess.CalledProcessError) as e:
            log.warning("Per-site dN/dS plot failed: %s", e)
    else:
        log.warning("%s missing — skipping dN/dS plot.", rst)

    # 3) Structure + PyMOL script
    pdb_path, _ = fetch_structure(target, outdir, prefix, pdb_id, uniprot=uniprot)
    if pdb_path is None:
        log.warning("No structure available — skipping PyMOL script.")
        return

    write_pymol_script(
        outdir / f"{prefix}_sites.pml",
        pdb_path,
        sites,
        cds_mapping,
        chain=chain,
        resi_offset=resi_offset,
    )


# =====================================================================
#  Pipeline
# =====================================================================


def symbol_to_prefix(symbol: str) -> str:
    """Derive a filesystem-friendly filename stem from a gene symbol.

    The :mod:`scripts.fetch_ensembl_msa_tree` filenames use underscores,
    so ``HLA-DQB1`` → ``HLA_DQB1`` and any other non-alphanumeric
    character is replaced with ``_``.

    Args:
        symbol: HGNC gene symbol (e.g. ``"HLA-DQB1"`` or ``"TP53"``).

    Returns:
        A sanitised stem suitable for use as ``--prefix``.
    """
    return re.sub(r"[^A-Za-z0-9_]+", "_", symbol).strip("_")


def run_pipeline(
    gene_symbol: str,
    outdir: Path | None,
    taxon: int,
    skip_codeml: bool,
    skip_post_analysis: bool = False,
    uniprot: str | None = None,
    pdb_id: str | None = None,
    chain: str = "A",
    resi_offset: int = 0,
    plot_threshold: float = 0.50,
) -> None:
    """Prepare inputs, write CTL files, run CodeML and report LRTs.

    Resolves the gene symbol (e.g. ``HLA-DQB1``) to its Ensembl gene,
    canonical protein stable ID and UniProt accession via
    :func:`resolve_ids.resolve` (shared with the ``protein_characteriser``
    pipeline), then drives the CodeML site-model workflow.

    Args:
        gene_symbol: HGNC gene symbol (e.g. ``"HLA-DQB1"``).
        outdir: Working directory. When ``None``, defaults to
            ``Path(symbol_to_prefix(gene_symbol))``.
        taxon: NCBI taxon ID passed through to the fetch script.
        skip_codeml: If ``True``, write the CTL files but skip execution
            (useful for inspecting the generated controls).
        skip_post_analysis: If ``True``, stop after the M8a-M8 LRT and
            don't extract BEB sites, plot dN/dS or build the PyMOL script.
        uniprot: Explicit UniProt accession to use for the AlphaFold
            fetch and TSV metadata. Overrides the Ensembl cross-reference
            lookup, which can return an isoform-level TrEMBL entry.
        pdb_id: Optional PDB accession for :func:`run_post_analysis`;
            when ``None`` an AlphaFold model is fetched instead.
        chain: Chain identifier of the target in the loaded structure.
        resi_offset: Residue-numbering offset applied inside the PyMOL
            script (mirrors the tutorial's ``alter ... resi=int(resi)+N``).
        plot_threshold: BEB probability cutoff for the dN/dS plot.
    """
    # ── Resolve identifiers via protein_characteriser ───────────────
    log.info("Resolving %r via protein_characteriser resolver …", gene_symbol)
    ids = resolve_identifier(gene_symbol, "gene")
    gene = ids["ensembl_gene"]
    target = ids["ensembl_protein"]
    if not gene or not target:
        raise RuntimeError(
            f"Could not resolve Ensembl gene/protein for symbol {gene_symbol!r} "
            f"(got gene={gene!r}, protein={target!r})."
        )
    # A user-supplied --uniprot wins over the cross-reference lookup
    # (e.g. P01920 vs. Q5Y7D6 for HLA-DQB1).
    uniprot = uniprot or ids["uniprot"]
    prefix = symbol_to_prefix(ids["gene_symbol"] or gene_symbol)

    workdir: Path = outdir if outdir is not None else Path(prefix)

    log.info(
        "Resolved %s → gene=%s protein=%s uniprot=%s prefix=%s outdir=%s",
        gene_symbol,
        gene,
        target,
        uniprot,
        prefix,
        workdir,
    )

    workdir.mkdir(parents=True, exist_ok=True)

    phy = workdir / f"{prefix}_subset.cds.mafft.trimal.phy"
    tree = workdir / f"{prefix}_subset.tree"

    ensure_inputs(gene, prefix, target, workdir, taxon, phy, tree)

    ctl_m8 = workdir / f"{prefix}_M8.ctl"
    ctl_m8a = workdir / f"{prefix}_M8a.ctl"
    mlc_m8 = f"{prefix}_M8.mlc"
    mlc_m8a = f"{prefix}_M8a.mlc"

    # CTL files use paths relative to outdir (since codeml runs with cwd=outdir)
    write_ctl(
        ctl_m8,
        seqfile=phy.name,
        treefile=tree.name,
        outfile=mlc_m8,
        nssites="8",
        fix_omega=0,
    )
    write_ctl(
        ctl_m8a,
        seqfile=phy.name,
        treefile=tree.name,
        outfile=mlc_m8a,
        nssites="8",
        fix_omega=1,
    )

    if skip_codeml:
        log.info("--skip-codeml set – stopping after CTL generation.")
        return

    # Run M8 and M8a in parallel — CodeML is single-processor and each
    # run gets its own subdirectory (see run_codeml), so they don't
    # interfere with each other's scratch files.
    mlc_m8_path = workdir / mlc_m8
    mlc_m8a_path = workdir / mlc_m8a
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        fut_m8 = pool.submit(run_codeml, ctl_m8, workdir, mlc=mlc_m8_path)
        fut_m8a = pool.submit(run_codeml, ctl_m8a, workdir, mlc=mlc_m8a_path)
        fut_m8.result()
        fut_m8a.result()

    # Preserve the M8 rst for the BEB/site-level analysis.  Each model
    # ran in its own subdirectory and scratch files were prefixed with
    # the CTL stem on move-back, so M8's rst is ``{prefix}_M8_rst``.
    rst_src = workdir / f"{prefix}_M8_rst"
    rst_dst = workdir / f"{prefix}_M8.rst.txt"
    if rst_src.exists() and (
        not rst_dst.exists() or rst_src.stat().st_mtime >= rst_dst.stat().st_mtime
    ):
        shutil.move(str(rst_src), str(rst_dst))
        log.info("Moved %s → %s", rst_src.name, rst_dst.name)
    elif not rst_dst.exists():
        log.warning("No rst file produced by CodeML – BEB analysis will be unavailable.")

    report_lrts(prefix, workdir)

    if skip_post_analysis:
        log.info("--skip-post-analysis set – stopping after LRT.")
        return

    run_post_analysis(
        prefix=prefix,
        target=target,
        outdir=workdir,
        uniprot=uniprot,
        pdb_id=pdb_id,
        chain=chain,
        resi_offset=resi_offset,
        plot_threshold=plot_threshold,
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    ap = argparse.ArgumentParser(
        description="Run CodeML site-models M8 and M8a for pervasive "
        "positive selection, following the evosite3d tutorial.",
    )
    ap.add_argument(
        "--gene-symbol",
        help="HGNC gene symbol (e.g. 'HLA-DQB1', 'TP53'). Resolved to "
        "Ensembl gene, canonical protein and UniProt via the shared "
        "software/protein_characteriser/bin/resolve_ids.py helper.",
    )
    ap.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Working directory for inputs, CTL files and CodeML "
        "outputs. Defaults to a sanitised form of the gene "
        "symbol (e.g. HLA-DQB1 → ./HLA_DQB1).",
    )
    ap.add_argument(
        "--taxon",
        type=int,
        default=9347,
        help="NCBI taxon ID for the Compara gene-tree prune "
        "(default 9347 = Boreoeutheria). Pass 7742 "
        "(Vertebrata) to make the full reference-species "
        "panel spanning fishes to mammals available to the "
        "species filter.",
    )
    ap.add_argument(
        "--skip-codeml", action="store_true", help="Only generate CTL files; do not run codeml."
    )
    ap.add_argument(
        "--skip-post-analysis",
        action="store_true",
        help="Stop after the M8a-M8 LRT: don't extract BEB "
        "sites, plot dN/dS or generate a PyMOL script.",
    )
    ap.add_argument(
        "--uniprot",
        default=None,
        help="Override the Ensembl→UniProt cross-reference and "
        "use this accession for the AlphaFold fetch and "
        "TSV metadata (e.g. P01920 for HLA-DQB1, where the "
        "automatic lookup may resolve to an isoform-level "
        "TrEMBL entry).",
    )
    ap.add_argument(
        "--pdb",
        default=None,
        help="Override AlphaFold lookup and use this RCSB PDB "
        "accession for the PyMOL visualisation.",
    )
    ap.add_argument(
        "--chain",
        default="A",
        help="Chain ID of the target protein in the loaded "
        "structure (default A, matches AlphaFold).",
    )
    ap.add_argument(
        "--resi-offset",
        type=int,
        default=0,
        help="Offset added to CDS positions inside the PyMOL "
        "script to align with crystal numbering "
        "(e.g. 32 for the tutorial's 1UVQ example).",
    )
    ap.add_argument(
        "--plot-threshold",
        type=float,
        default=0.50,
        help="BEB probability threshold for the per-site dN/dS plot (default 0.50).",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    run_pipeline(
        gene_symbol=args.gene_symbol,
        outdir=args.outdir,
        taxon=args.taxon,
        skip_codeml=args.skip_codeml,
        skip_post_analysis=args.skip_post_analysis,
        uniprot=args.uniprot,
        pdb_id=args.pdb,
        chain=args.chain,
        resi_offset=args.resi_offset,
        plot_threshold=args.plot_threshold,
    )


if __name__ == "__main__":
    main()
