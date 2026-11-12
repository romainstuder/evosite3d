"""Shared utilities for the positive-selection site-model pipeline.

This module is imported by the three pipeline scripts:

    1. ``prepare_data.py``          — data acquisition and preparation
    2. ``run_codeml_site_models.py`` — CodeML execution and result extraction
    3. ``analyse_site_models.py``    — visualisation (Jalview, PyMOL, plots)

It provides path constants, ID resolution helpers, data classes for
parsed CodeML outputs, and common CLI argument definitions.
"""

from __future__ import annotations

import argparse
import logging
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ── Repository layout ─────────────────────────────────────────────────
_REPO = Path(__file__).resolve().parent.parent.parent
GET_POS_SCRIPT = _REPO / "scripts" / "get_position_cds_trimal.py"
PLOT_SCRIPT = _REPO / "scripts" / "plot_dnds_per_site.py"

# Make the shared resolver importable.
sys.path.insert(0, str(_REPO / "bin"))
sys.path.insert(0, str(_REPO / "lib"))
sys.path.insert(0, str(_REPO / "scripts"))

from resolve_ids import resolve as resolve_identifier  # noqa: E402  # type: ignore

# ── URLs ──────────────────────────────────────────────────────────────
ALPHAFOLD_FILES = "https://alphafold.ebi.ac.uk/files"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"
RCSB_FILES = "https://files.rcsb.org/download"

# ── Logging ───────────────────────────────────────────────────────────
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("positive_selection_detector")


# =====================================================================
#  Small helpers
# =====================================================================


def check_tool(name: str) -> None:
    """Raise ``RuntimeError`` if an external binary is missing from ``PATH``."""
    if shutil.which(name) is None:
        raise RuntimeError(f"Required tool '{name}' not found on PATH.")


def symbol_to_prefix(symbol: str) -> str:
    """Derive a filesystem-friendly filename stem from a gene symbol.

    ``HLA-DQB1`` → ``HLA_DQB1``; any non-alphanumeric character is
    replaced with ``_``.
    """
    return re.sub(r"[^A-Za-z0-9_]+", "_", symbol).strip("_")


def http_get(url: str, timeout: int = 60) -> bytes | None:
    """Simple HTTP GET returning raw bytes, or ``None`` on failure."""
    try:
        req = Request(url, headers={"User-Agent": "evosite3d/positive_selection_detector"})
        with urlopen(req, timeout=timeout) as resp:
            return resp.read()
    except (HTTPError, URLError, TimeoutError) as e:
        log.warning("HTTP GET failed for %s: %s", url, e)
        return None


# =====================================================================
#  FASTA helpers
# =====================================================================


def read_fasta_dict(in_fasta: Path) -> dict[str, str]:
    """Read a FASTA file into a ``{id: sequence}`` mapping."""
    return {r.id: str(r.seq) for r in SeqIO.parse(in_fasta, "fasta")}


def write_seq_fasta(seqs: dict[str, str], out_fasta: Path) -> None:
    """Write a ``{id: sequence}`` mapping as multi-FASTA."""
    records = [SeqRecord(Seq(seq), id=sid, description="") for sid, seq in seqs.items()]
    SeqIO.write(records, out_fasta, "fasta")


# =====================================================================
#  Data classes — CodeML parsed outputs
# =====================================================================


@dataclass
class LnL:
    """A single ``lnL`` line extracted from an ``.mlc`` file."""

    ntime: int
    np: int
    lnl: float


_LNL_RE = re.compile(r"lnL\(ntime:\s*(\d+)\s+np:\s*(\d+)\):\s*(-?\d+\.\d+)")


def parse_lnl(mlc: Path) -> list[LnL]:
    """Extract every ``lnL`` line from a CodeML ``.mlc`` file."""
    entries = []
    for line in mlc.read_text().splitlines():
        m = _LNL_RE.search(line)
        if m:
            entries.append(LnL(int(m.group(1)), int(m.group(2)), float(m.group(3))))
    return entries


@dataclass
class BebSite:
    """A positively-selected site from a CodeML BEB table."""

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
    """Extract the BEB ``Positively selected sites`` block from an ``.mlc``.

    Scans for the *last* ``Bayes Empirical Bayes (BEB) analysis`` header
    and reads ``<pos> <aa> <prob>[*|**] <mean> +- <se>`` rows until a
    blank line.
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
#  Map trimal column indices → target-protein positions
# =====================================================================


_GET_POS_RE = re.compile(r"site_index:\s*(\d+)\s+CDS pos:\s*(\d+|N/A)\s+AA:\s*(\S+)")


def map_sites_to_cds(
    full_aln: Path,
    trimal_cols: Path,
    target: str,
    sites: list[int],
) -> dict[int, tuple[int | None, str]]:
    """Map trimal column indices to CDS positions in ``target``.

    Shells out to ``scripts/get_position_cds_trimal.py`` and parses
    its ``site_index:`` lines. Returns ``{trimal_pos: (cds_pos, aa)}``.
    """
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


# =====================================================================
#  ID resolution + workdir setup
# =====================================================================


@dataclass
class ResolvedIDs:
    """Resolved identifiers and paths for a gene symbol."""

    gene: str
    target: str
    uniprot: str | None
    prefix: str
    workdir: Path


def resolve_gene(
    gene_symbol: str,
    outdir: Path | None = None,
    uniprot_override: str | None = None,
) -> ResolvedIDs:
    """Resolve a gene symbol to Ensembl/UniProt IDs and working directory.

    Args:
        gene_symbol: HGNC gene symbol (e.g. ``"HLA-DQB1"``).
        outdir: Working directory override.  Defaults to ``./PREFIX``.
        uniprot_override: Explicit UniProt accession (overrides lookup).

    Returns:
        A :class:`ResolvedIDs` with all fields populated.
    """
    log.info("Resolving %r via protein_characteriser resolver …", gene_symbol)
    ids = resolve_identifier(gene_symbol, "gene")
    gene = ids["ensembl_gene"]
    target = ids["ensembl_protein"]
    if not gene or not target:
        raise RuntimeError(
            f"Could not resolve Ensembl gene/protein for symbol {gene_symbol!r} "
            f"(got gene={gene!r}, protein={target!r})."
        )
    uniprot = uniprot_override or ids["uniprot"]
    prefix = symbol_to_prefix(ids["gene_symbol"] or gene_symbol)
    workdir = outdir if outdir is not None else Path(prefix)

    log.info(
        "Resolved %s → gene=%s protein=%s uniprot=%s prefix=%s outdir=%s",
        gene_symbol,
        gene,
        target,
        uniprot,
        prefix,
        workdir,
    )
    return ResolvedIDs(gene=gene, target=target, uniprot=uniprot, prefix=prefix, workdir=workdir)


# =====================================================================
#  Common CLI arguments
# =====================================================================


def add_common_args(ap: argparse.ArgumentParser) -> None:
    """Add arguments shared by all three pipeline scripts."""
    ap.add_argument(
        "--gene-symbol",
        required=True,
        help="HGNC gene symbol (e.g. 'HLA-DQB1', 'TP53').",
    )
    ap.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Working directory. Defaults to a sanitised form of the gene symbol.",
    )
    ap.add_argument(
        "--uniprot",
        default=None,
        help="Override the Ensembl→UniProt cross-reference.",
    )
