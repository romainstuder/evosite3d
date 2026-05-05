#!/usr/bin/env python3
"""
modules/pervasive_selection.py
===============================
Per-residue pervasive positive selection analysis using PAML/CodeML
site models.

  +--------------------------------------------------------------+
  |  PERVASIVE POSITIVE SELECTION (CodeML Site Models)           |
  |                                                              |
  |  1. Fetch CDS sequences from Ensembl for orthologs in MSA   |
  |  2. Build codon alignment (protein MSA -> back-translate)    |
  |  3. Trim with trimal, convert to PHYLIP                     |
  |  4. Fetch + prune gene tree to match alignment taxa          |
  |  5. Run CodeML site models:                                  |
  |     M0, M1a, M2a, M3, M7, M8, M8a                          |
  |  6. Likelihood Ratio Tests:                                  |
  |     - M0 vs M3  (heterogeneous rates)                       |
  |     - M1a vs M2a (positive selection, conservative)         |
  |     - M7 vs M8  (positive selection, beta model)            |
  |     - M8 vs M8a (positive selection, preferred test)        |
  |  7. Parse BEB (Bayes Empirical Bayes) for sites under       |
  |     positive selection from M2a and M8                       |
  |  8. Map trimmed alignment positions -> human sequence        |
  |                                                              |
  |  Falls back to empty results when codeml is unavailable.    |
  +--------------------------------------------------------------+

Inputs:  resolved_ids.json, sequence.fasta, orthologs.fasta (protein MSA)
Outputs: pervasive_selection.tsv
         codon_alignment.fasta (intermediate)
         codeml_results.json   (LRT summary)
"""

import argparse
import csv
import json
import logging
import math
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from lib.ensembl import extract_protein_ids_from_msa, fetch_cds_sequences, fetch_gene_tree_newick

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


# =====================================================================
#  2. Build codon alignment
# =====================================================================


def _build_codon_alignment(
    protein_msa_path: Path,
    cds_seqs: dict[str, str],
) -> MultipleSeqAlignment | None:
    """Create a codon alignment by back-translating a protein MSA.

    For each sequence in the protein MSA that has a matching CDS, maps
    each aligned amino acid to its corresponding codon and gaps to ``---``.

    Args:
        protein_msa_path: Path to the aligned protein FASTA.
        cds_seqs: Dict mapping sequence IDs to CDS nucleotide strings.

    Returns:
        A :class:`MultipleSeqAlignment` of codon-aligned nucleotides,
        or ``None`` if fewer than 4 sequences could be aligned.
    """
    protein_aln = AlignIO.read(protein_msa_path, "fasta")
    aligned_records = []

    for prot_rec in protein_aln:
        # Find matching CDS - try exact id, then search for ensembl protein ID
        nt_seq = None
        rec_id = prot_rec.id

        if rec_id in cds_seqs:
            nt_seq = cds_seqs[rec_id]
        else:
            match = re.search(r"(ENS\w*P\d+)", rec_id)
            if match:
                ens_id = match.group(1)
                if ens_id in cds_seqs:
                    nt_seq = cds_seqs[ens_id]
                    rec_id = ens_id

        if nt_seq is None:
            continue

        prot_aligned = str(prot_rec.seq)
        prot_nogaps = prot_aligned.replace("-", "")

        # Check length compatibility (CDS should be ~3x protein, +/- stop codon)
        expected_len = len(prot_nogaps) * 3
        if not (expected_len <= len(nt_seq) <= expected_len + 3):
            log.debug(
                "Length mismatch for %s: protein=%d AA, CDS=%d nt (expected %d)",
                rec_id,
                len(prot_nogaps),
                len(nt_seq),
                expected_len,
            )
            continue

        # Back-translate: map protein alignment gaps to codon gaps
        codon_aligned = []
        nt_pos = 0
        for aa in prot_aligned:
            if aa == "-":
                codon_aligned.append("---")
            else:
                if nt_pos + 3 <= len(nt_seq):
                    codon_aligned.append(nt_seq[nt_pos : nt_pos + 3])
                else:
                    codon_aligned.append("NNN")
                nt_pos += 3

        aligned_records.append(
            SeqRecord(
                Seq("".join(codon_aligned)),
                id=prot_rec.id,
                description="",
            )
        )

    if len(aligned_records) < 4:
        log.warning("Only %d sequences with CDS – need at least 4.", len(aligned_records))
        return None

    log.info("Codon alignment: %d sequences.", len(aligned_records))
    return MultipleSeqAlignment(aligned_records)


# =====================================================================
#  3. Trim alignment and convert to PHYLIP
# =====================================================================


def _trim_alignment(
    codon_aln: MultipleSeqAlignment,
    workdir: Path,
) -> tuple[Path | None, list[int] | None]:
    """Trim a codon alignment with trimal and return column mapping.

    Runs ``trimal -automated1`` if available, otherwise returns the
    untrimmed alignment. Also produces a column numbering file for
    mapping trimmed positions back to the original alignment.

    Args:
        codon_aln: The codon :class:`MultipleSeqAlignment`.
        workdir: Temporary working directory for intermediate files.

    Returns:
        Tuple of (path to trimmed FASTA, list of kept original column
        indices) or ``(None, None)`` on failure.
    """
    input_fasta = workdir / "codon_aln.fasta"
    AlignIO.write(codon_aln, input_fasta, "fasta")

    trimal = shutil.which("trimal")
    if not trimal:
        log.info("trimal not found – using untrimmed alignment.")
        # Return identity column mapping (codon-based: every 3 nt = 1 codon)
        aln_len = codon_aln.get_alignment_length()
        cols = list(range(0, aln_len, 3))
        return input_fasta, cols

    trimmed_fasta = workdir / "codon_aln_trimmed.fasta"
    cols_file = workdir / "trimmed_cols.txt"

    cmd = [
        trimal,
        "-in",
        str(input_fasta),
        "-out",
        str(trimmed_fasta),
        "-gt",
        "0.9",
        "-resoverlap",
        "0.75",
        "-seqoverlap",
        "85",
        "-colnumbering",
        "-block",
        "3",
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
        )
        # trimal writes column numbering to stdout
        if result.stdout.strip():
            cols_file.write_text(result.stdout.strip())
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        log.warning("trimal failed: %s – using untrimmed.", e)
        aln_len = codon_aln.get_alignment_length()
        cols = list(range(0, aln_len, 3))
        return input_fasta, cols

    if not trimmed_fasta.exists() or trimmed_fasta.stat().st_size == 0:
        log.warning("trimal produced no output – using untrimmed.")
        aln_len = codon_aln.get_alignment_length()
        cols = list(range(0, aln_len, 3))
        return input_fasta, cols

    # Parse kept columns from trimal output
    kept_cols = None
    if cols_file.exists():
        text = cols_file.read_text().strip()
        # trimal outputs comma-separated column indices
        try:
            kept_cols = [int(c.strip()) for c in text.split(",") if c.strip()]
        except ValueError:
            kept_cols = None

    if kept_cols is None:
        # Fallback: assume all columns kept
        try:
            trimmed = AlignIO.read(trimmed_fasta, "fasta")
            kept_cols = list(range(trimmed.get_alignment_length()))
        except Exception:
            kept_cols = []

    log.info("Trimmed alignment: %d columns kept.", len(kept_cols))
    return trimmed_fasta, kept_cols


def _fasta_to_phylip(fasta_path: Path, phylip_path: Path, name_length: int = 50):
    """Convert a FASTA alignment to PHYLIP sequential format.

    Produces PHYLIP output compatible with CodeML, using sequence names
    truncated to ``name_length`` characters.

    Args:
        fasta_path: Input aligned FASTA file.
        phylip_path: Output PHYLIP file path.
        name_length: Maximum length for sequence names (default 50).
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return

    n_seqs = len(records)
    aln_len = len(records[0].seq)

    with open(phylip_path, "w") as f:
        f.write(f"{n_seqs}\t{aln_len}\n")
        for rec in records:
            name = rec.id[:name_length].rstrip("_").replace(" ", "")
            f.write(f"{name}  {str(rec.seq)}\n")


# =====================================================================
#  4. Fetch and prune gene tree
# =====================================================================


def _prune_tree_to_taxa(newick: str, taxa: list[str], workdir: Path) -> Path | None:
    """Prune a Newick tree to keep only the specified taxa.

    Uses ``nw_prune -v`` (Newick Utilities) if available. Falls back to
    a simple regex-based approach that extracts matching leaf labels.

    The output tree file includes the ``n_taxa 1`` header line required
    by CodeML.

    Args:
        newick: Full Newick tree string.
        taxa: List of taxon/leaf names to retain.
        workdir: Temporary directory for output file.

    Returns:
        Path to the pruned tree file, or ``None`` on failure.
    """
    tree_path = workdir / "pruned.tree"

    nw_prune = shutil.which("nw_prune")
    if nw_prune:
        try:
            # nw_prune -v keeps only listed taxa (inverted prune)
            cmd = [nw_prune, "-v", "-", *taxa]
            result = subprocess.run(
                cmd,
                input=newick,
                capture_output=True,
                text=True,
                timeout=30,
            )
            if result.returncode == 0 and result.stdout.strip():
                with open(tree_path, "w") as f:
                    f.write(f"{len(taxa)} 1\n")
                    f.write(result.stdout.strip() + "\n")
                return tree_path
        except Exception as e:
            log.warning("nw_prune failed: %s", e)

    # Fallback: write the full tree with taxa count header
    # CodeML can sometimes handle extra taxa, but results may be suboptimal
    log.info("nw_prune not available – writing full tree (CodeML may warn).")
    with open(tree_path, "w") as f:
        f.write(f"{len(taxa)} 1\n")
        f.write(newick + "\n")
    return tree_path


# =====================================================================
#  5. Run CodeML
# =====================================================================

_CTL_TEMPLATE = """\
     seqfile = {seqfile}
    treefile = {treefile}
     outfile = {outfile}

       noisy = 0
     verbose = 1
     runmode = 0

     seqtype = 1
   CodonFreq = 2
       clock = 0
      aaDist = 0
       model = 0
     NSsites = {nssites}
       icode = 0
       Mgene = 0
   fix_kappa = 0
       kappa = 2
   fix_omega = {fix_omega}
       omega = {omega}
       getSE = 0
RateAncestor = 0
  Small_Diff = .45e-6
   cleandata = 0
 fix_blength = 0
"""


def _run_codeml(
    phylip_path: Path,
    tree_path: Path,
    workdir: Path,
    nssites: str = "0 1 2 3 7 8",
    fix_omega: int = 0,
    omega: float = 1.0,
    label: str = "main",
) -> tuple[Path | None, Path | None]:
    """Run CodeML with the given site models.

    Writes a control file, executes ``codeml``, and returns paths to the
    main output (``.mlc``) and the ``rst`` file containing BEB results.

    Args:
        phylip_path: Path to the PHYLIP alignment.
        tree_path: Path to the pruned Newick tree.
        workdir: Working directory for CodeML execution.
        nssites: Space-separated NSsites model numbers.
        fix_omega: Whether to fix omega (1) or estimate it (0).
        omega: Initial or fixed omega value.
        label: Label for output file naming.

    Returns:
        Tuple of (mlc_path, rst_path) or ``(None, None)`` on failure.
    """
    codeml = shutil.which("codeml")
    if not codeml:
        log.warning("codeml not found in PATH – cannot run site models.")
        return None, None

    mlc_path = workdir / f"{label}.mlc"
    ctl_path = workdir / f"{label}.ctl"
    rst_path = workdir / "rst"

    ctl_content = _CTL_TEMPLATE.format(
        seqfile=phylip_path.name,
        treefile=tree_path.name,
        outfile=mlc_path.name,
        nssites=nssites,
        fix_omega=fix_omega,
        omega=omega,
    )
    ctl_path.write_text(ctl_content)

    log.info("Running codeml (NSsites=%s, label=%s) ...", nssites, label)
    try:
        result = subprocess.run(
            [codeml, ctl_path.name],
            cwd=str(workdir),
            capture_output=True,
            text=True,
            timeout=5400,  # 90 min
        )
        if result.returncode != 0:
            log.warning("codeml exited with code %d: %s", result.returncode, result.stderr[:500])
    except subprocess.TimeoutExpired:
        log.warning("codeml timed out after 90 minutes.")
        return None, None
    except FileNotFoundError:
        log.warning("codeml not found.")
        return None, None

    if not mlc_path.exists():
        log.warning("codeml did not produce %s", mlc_path)
        return None, None

    rst_out = workdir / f"{label}.rst"
    if rst_path.exists():
        shutil.copy2(rst_path, rst_out)

    return mlc_path, rst_out if rst_out.exists() else None


# =====================================================================
#  6. Parse CodeML output
# =====================================================================


def _parse_lnl_values(mlc_path: Path) -> list[dict]:
    """Parse log-likelihood values from a CodeML mlc file.

    Extracts ``lnL``, number of parameters (``np``), and the model
    number from each ``lnL(ntime: ... np: ...):`` line.

    Args:
        mlc_path: Path to the ``.mlc`` output file.

    Returns:
        List of dicts with keys ``np`` and ``lnL``, in the order they
        appear in the file (matching the NSsites order).
    """
    results = []
    text = mlc_path.read_text()

    for m in re.finditer(r"lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s+([-\d.]+)", text):
        np_val = int(m.group(1))
        lnl_val = float(m.group(2))
        results.append({"np": np_val, "lnL": lnl_val})

    return results


def _compute_lrt(lnl_null: float, lnl_alt: float, df: int) -> dict:
    """Compute a Likelihood Ratio Test statistic and p-value.

    Uses the chi-squared approximation: ``2 * (lnL_alt - lnL_null)``
    with ``df`` degrees of freedom.

    Args:
        lnl_null: Log-likelihood of the null model.
        lnl_alt: Log-likelihood of the alternative model.
        df: Degrees of freedom (difference in number of parameters).

    Returns:
        Dict with ``lrt_statistic``, ``df``, ``p_value``, and
        ``significant`` (at alpha=0.05).
    """
    lrt = 2 * (lnl_alt - lnl_null)
    lrt = max(0, lrt)  # Cannot be negative

    # Chi-squared survival function approximation
    # Using scipy if available, else a simple approximation
    try:
        from scipy.stats import chi2

        p_value = chi2.sf(lrt, df)
    except ImportError:
        # Wilson-Hilferty approximation for chi2 CDF
        p_value = _chi2_sf_approx(lrt, df)

    return {
        "lrt_statistic": round(lrt, 4),
        "df": df,
        "p_value": p_value,
        "significant": p_value < 0.05,
    }


def _chi2_sf_approx(x: float, df: int) -> float:
    """Approximate chi-squared survival function without scipy.

    Uses a normal approximation for large degrees of freedom and a
    direct series for small df values.

    Args:
        x: Chi-squared test statistic.
        df: Degrees of freedom.

    Returns:
        Approximate p-value (survival function value).
    """
    if x <= 0:
        return 1.0
    if df <= 0:
        return 0.0

    # For df=1 and df=2, use exact formulas
    if df == 1:
        z = math.sqrt(x)
        return math.erfc(z / math.sqrt(2))
    if df == 2:
        return math.exp(-x / 2)

    # Wilson-Hilferty normal approximation
    z = ((x / df) ** (1 / 3) - (1 - 2 / (9 * df))) / math.sqrt(2 / (9 * df))
    # Standard normal survival function
    return 0.5 * math.erfc(z / math.sqrt(2))


def _parse_beb_from_mlc(mlc_path: Path) -> list[dict]:
    """Parse Bayes Empirical Bayes (BEB) results from a CodeML mlc file.

    Looks for the BEB section under M2a and M8 models and extracts
    per-site posterior probabilities of positive selection.

    Args:
        mlc_path: Path to the ``.mlc`` output file.

    Returns:
        List of dicts with keys ``aln_position`` (1-based position in
        the trimmed alignment), ``aa``, ``beb_prob`` (Pr(w>1)),
        ``post_mean_omega``, ``post_se_omega``, and ``model``.
    """
    beb_sites = []
    text = mlc_path.read_text()

    # Find BEB sections
    beb_pattern = re.compile(
        r"Bayes Empirical Bayes.*?Positively selected sites.*?\n"
        r"\(amino acids refer to 1st sequence:.*?\)\n\n"
        r"\s+Pr\(w>1\).*?\n\n(.*?)(?:\n\n|\Z)",
        re.DOTALL,
    )

    for beb_match in beb_pattern.finditer(text):
        block = beb_match.group(1)
        for line in block.strip().split("\n"):
            line = line.strip()
            if not line:
                continue
            # Format: "  8 R      0.971*        2.453 +- 0.275"
            m = re.match(
                r"\s*(\d+)\s+([A-Z*])\s+([\d.]+)\*{0,2}\s+([\d.]+)\s*\+-\s*([\d.]+)",
                line,
            )
            if m:
                beb_sites.append(
                    {
                        "aln_position": int(m.group(1)),
                        "aa": m.group(2),
                        "beb_prob": float(m.group(3)),
                        "post_mean_omega": float(m.group(4)),
                        "post_se_omega": float(m.group(5)),
                    }
                )

    return beb_sites


def _parse_beb_from_rst(rst_path: Path) -> list[dict]:
    """Parse per-site dN/dS and class probabilities from the rst file.

    The rst file contains detailed per-site BEB output at the end,
    with columns: position, AA, class probabilities, most likely class,
    estimated dN/dS +/- SE.

    Args:
        rst_path: Path to the CodeML ``rst`` output file.

    Returns:
        List of dicts with ``aln_position``, ``aa``,
        ``estimated_omega``, ``omega_se``, and ``most_likely_class``.
    """
    sites = []
    if not rst_path or not rst_path.exists():
        return sites

    text = rst_path.read_text()

    # Find per-site BEB section at end of rst
    # Format: "   1 M   0.01184 ... ( 3)  0.352 +-  0.154"
    pattern = re.compile(
        r"^\s*(\d+)\s+([A-Z*])\s+[\d.]+(?:\s+[\d.]+)*\s+\(\s*(\d+)\)\s+([\d.]+)\s*\+-\s*([\d.]+)",
        re.MULTILINE,
    )

    for m in pattern.finditer(text):
        sites.append(
            {
                "aln_position": int(m.group(1)),
                "aa": m.group(2),
                "most_likely_class": int(m.group(3)),
                "estimated_omega": float(m.group(4)),
                "omega_se": float(m.group(5)),
            }
        )

    return sites


# =====================================================================
#  7. Position mapping
# =====================================================================


def _build_aln_to_human_map(
    msa_path: Path,
    human_id: str | None,
) -> dict[int, int]:
    """Map trimmed alignment codon positions to human sequence positions.

    Walks through the human sequence in the MSA, counting non-gap
    positions to build a mapping from alignment column (1-based codon
    index) to human sequence position (1-based amino acid index).

    Args:
        msa_path: Path to the codon alignment FASTA.
        human_id: Identifier of the human sequence in the MSA.

    Returns:
        Dict mapping alignment position (1-based) to human AA position
        (1-based).
    """
    mapping = {}
    records = list(SeqIO.parse(msa_path, "fasta"))
    if not records:
        return mapping

    # Find human sequence
    human_rec = None
    if human_id:
        for rec in records:
            if human_id in rec.id or rec.id in human_id:
                human_rec = rec
                break
    if not human_rec:
        for rec in records:
            if "homo_sapiens" in rec.id.lower() or "human" in rec.id.lower():
                human_rec = rec
                break
    if not human_rec:
        human_rec = records[0]  # fallback: first sequence

    seq = str(human_rec.seq)
    aln_len = len(seq)

    # For codon alignment, step by 3
    human_pos = 0
    codon_idx = 0
    for i in range(0, aln_len, 3):
        codon_idx += 1
        codon = seq[i : i + 3]
        if codon != "---" and "N" not in codon:
            human_pos += 1
            mapping[codon_idx] = human_pos

    return mapping


# =====================================================================
#  PUBLIC PIPELINE
# =====================================================================

COLUMNS = [
    "position",
    "aa",
    "codeml_beb_prob",
    "codeml_post_mean_omega",
    "codeml_post_se_omega",
    "codeml_estimated_omega",
    "codeml_omega_se",
    "codeml_selection_class",
    "codeml_lrt_m1a_m2a_pval",
    "codeml_lrt_m7_m8_pval",
    "codeml_lrt_m8_m8a_pval",
    "codeml_significant",
]


def annotate_pervasive_selection(
    ids: dict,
    sequence: str,
    msa_path: Path,
    outdir: Path,
) -> list[dict]:
    """Run the full pervasive selection analysis pipeline.

    Orchestrates CDS fetching, codon alignment building, tree pruning,
    CodeML execution, and result parsing.

    Args:
        ids: Resolved identifiers dict (from ``resolved_ids.json``).
        sequence: Human protein sequence string.
        msa_path: Path to the protein MSA FASTA from the evolution module.
        outdir: Output directory for intermediate and result files.

    Returns:
        List of per-residue dicts with CodeML site model results.
    """
    seq_len = len(sequence)
    results = [{"position": i + 1, "aa": sequence[i]} for i in range(seq_len)]

    ensembl_gene = ids.get("ensembl_gene", "")
    if not ensembl_gene:
        log.warning("No Ensembl gene ID – skipping pervasive selection.")
        return results

    if not msa_path.exists() or msa_path.stat().st_size == 0:
        log.warning("No MSA available – skipping pervasive selection.")
        return results

    # Check codeml availability early
    codeml = shutil.which("codeml")
    if not codeml:
        log.warning(
            "codeml (PAML) not found in PATH – skipping pervasive selection. "
            "Install PAML: http://abacus.gene.ucl.ac.uk/software/paml.html"
        )
        return results

    # ── Step 1: Fetch CDS sequences ───────────────────────────────
    protein_ids = extract_protein_ids_from_msa(msa_path)
    if len(protein_ids) < 4:
        log.warning("Only %d Ensembl protein IDs in MSA – need >= 4.", len(protein_ids))
        return results

    cds_seqs = fetch_cds_sequences(ensembl_gene, protein_ids)
    if len(cds_seqs) < 4:
        log.warning("Only %d CDS sequences fetched – need >= 4.", len(cds_seqs))
        return results

    # ── Step 2: Build codon alignment ─────────────────────────────
    codon_aln = _build_codon_alignment(msa_path, cds_seqs)
    if codon_aln is None:
        return results

    # Save codon alignment
    codon_aln_path = outdir / "codon_alignment.fasta"
    AlignIO.write(codon_aln, codon_aln_path, "fasta")
    log.info("Wrote codon alignment: %s", codon_aln_path)

    # ── Step 3: Trim + convert to PHYLIP ──────────────────────────
    with tempfile.TemporaryDirectory(prefix="codeml_") as tmpdir:
        workdir = Path(tmpdir)

        trimmed_fasta, kept_cols = _trim_alignment(codon_aln, workdir)
        if trimmed_fasta is None:
            return results

        phylip_path = workdir / "alignment.phy"
        _fasta_to_phylip(trimmed_fasta, phylip_path)

        if not phylip_path.exists() or phylip_path.stat().st_size == 0:
            log.warning("PHYLIP conversion failed.")
            return results

        # ── Step 4: Fetch + prune tree ────────────────────────────
        newick = fetch_gene_tree_newick(ensembl_gene)
        if not newick:
            log.warning("Could not fetch gene tree – cannot run CodeML.")
            return results

        # Get taxa names from trimmed alignment
        trimmed_records = list(SeqIO.parse(trimmed_fasta, "fasta"))
        taxa = [rec.id for rec in trimmed_records]

        tree_path = _prune_tree_to_taxa(newick, taxa, workdir)
        if tree_path is None:
            log.warning("Tree pruning failed.")
            return results

        # Move human sequence to top of PHYLIP file (CodeML uses 1st seq as ref)
        _move_human_to_top(phylip_path, ids.get("ensembl_protein", ""))

        # ── Step 5: Run CodeML ────────────────────────────────────
        # Run M0, M1a, M2a, M3, M7, M8
        mlc_main, rst_main = _run_codeml(
            phylip_path,
            tree_path,
            workdir,
            nssites="0 1 2 3 7 8",
            label="site_models",
        )

        # Run M8a separately (fix_omega=1, omega=1)
        mlc_m8a, rst_m8a = _run_codeml(
            phylip_path,
            tree_path,
            workdir,
            nssites="8",
            fix_omega=1,
            omega=1.0,
            label="m8a",
        )

        # ── Step 6: Parse results ─────────────────────────────────
        lrt_results = {}

        if mlc_main:
            lnl_values = _parse_lnl_values(mlc_main)
            # Order: M0, M1a, M2a, M3, M7, M8
            if len(lnl_values) >= 6:
                m0, m1a, m2a, m3, m7, m8 = lnl_values[:6]

                lrt_results["M0_M3"] = _compute_lrt(
                    m0["lnL"],
                    m3["lnL"],
                    m3["np"] - m0["np"],
                )
                lrt_results["M1a_M2a"] = _compute_lrt(
                    m1a["lnL"],
                    m2a["lnL"],
                    m2a["np"] - m1a["np"],
                )
                lrt_results["M7_M8"] = _compute_lrt(
                    m7["lnL"],
                    m8["lnL"],
                    m8["np"] - m7["np"],
                )

                # M8 vs M8a
                if mlc_m8a:
                    lnl_m8a = _parse_lnl_values(mlc_m8a)
                    if lnl_m8a:
                        lrt_results["M8_M8a"] = _compute_lrt(
                            lnl_m8a[0]["lnL"],
                            m8["lnL"],
                            m8["np"] - lnl_m8a[0]["np"],
                        )

                log.info(
                    "LRT results: %s",
                    {
                        k: f"p={v['p_value']:.2e} {'*' if v['significant'] else ''}"
                        for k, v in lrt_results.items()
                    },
                )
            else:
                log.warning("Expected 6 lnL values, got %d.", len(lnl_values))

        # Save LRT summary
        lrt_summary = {
            k: {kk: (str(vv) if isinstance(vv, float) else vv) for kk, vv in v.items()}
            for k, v in lrt_results.items()
        }
        lrt_json_path = outdir / "codeml_results.json"
        with open(lrt_json_path, "w") as f:
            json.dump(lrt_summary, f, indent=2, default=str)
        log.info("Wrote LRT summary: %s", lrt_json_path)

        # ── Step 7: Parse BEB sites ───────────────────────────────
        beb_sites = {}  # aln_position -> best BEB info

        if mlc_main:
            beb_from_mlc = _parse_beb_from_mlc(mlc_main)
            for site in beb_from_mlc:
                pos = site["aln_position"]
                if pos not in beb_sites or site["beb_prob"] > beb_sites[pos]["beb_prob"]:
                    beb_sites[pos] = site

        rst_sites = {}  # aln_position -> omega estimate from rst
        if rst_main:
            rst_data = _parse_beb_from_rst(rst_main)
            for site in rst_data:
                rst_sites[site["aln_position"]] = site

        # ── Step 8: Map to human positions ────────────────────────
        aln_to_human = _build_aln_to_human_map(trimmed_fasta, ids.get("ensembl_protein"))

        # Determine significance
        any_significant = any(
            v.get("significant", False)
            for k, v in lrt_results.items()
            if k in ("M1a_M2a", "M7_M8", "M8_M8a")
        )

        # Global LRT p-values (for all residues)
        m1a_m2a_pval = lrt_results.get("M1a_M2a", {}).get("p_value", "")
        m7_m8_pval = lrt_results.get("M7_M8", {}).get("p_value", "")
        m8_m8a_pval = lrt_results.get("M8_M8a", {}).get("p_value", "")

        for r in results:
            r["codeml_lrt_m1a_m2a_pval"] = m1a_m2a_pval
            r["codeml_lrt_m7_m8_pval"] = m7_m8_pval
            r["codeml_lrt_m8_m8a_pval"] = m8_m8a_pval
            r["codeml_significant"] = any_significant

        # Assign per-site BEB scores
        n_beb = 0
        for aln_pos, beb in beb_sites.items():
            human_pos = aln_to_human.get(aln_pos)
            if human_pos is None or human_pos > seq_len:
                continue
            r = results[human_pos - 1]
            r["codeml_beb_prob"] = round(beb["beb_prob"], 4)
            r["codeml_post_mean_omega"] = round(beb["post_mean_omega"], 4)
            r["codeml_post_se_omega"] = round(beb["post_se_omega"], 4)
            n_beb += 1

            # Classification
            if beb["beb_prob"] >= 0.99:
                r["codeml_selection_class"] = "positive_99"
            elif beb["beb_prob"] >= 0.95:
                r["codeml_selection_class"] = "positive_95"
            elif beb["beb_prob"] >= 0.50:
                r["codeml_selection_class"] = "positive_50"

        # Assign per-site omega from rst
        for aln_pos, rst in rst_sites.items():
            human_pos = aln_to_human.get(aln_pos)
            if human_pos is None or human_pos > seq_len:
                continue
            r = results[human_pos - 1]
            r["codeml_estimated_omega"] = round(rst["estimated_omega"], 4)
            r["codeml_omega_se"] = round(rst["omega_se"], 4)

        n_95 = sum(
            1 for r in results if r.get("codeml_selection_class", "").startswith("positive_9")
        )
        log.info(
            "BEB: %d sites mapped, %d with Pr(w>1)>=0.95, overall significant=%s.",
            n_beb,
            n_95,
            any_significant,
        )

    return results


def _move_human_to_top(phylip_path: Path, ensembl_protein: str):
    """Move the human sequence to the first position in a PHYLIP file.

    CodeML uses the first sequence as the reference for BEB amino acid
    reporting, so the human sequence should be listed first.

    Args:
        phylip_path: Path to the PHYLIP file to reorder in-place.
        ensembl_protein: Ensembl protein ID to identify the human entry.
    """
    lines = phylip_path.read_text().strip().split("\n")
    if len(lines) < 3:
        return

    header = lines[0]
    seq_lines = lines[1:]

    human_idx = None
    for i, line in enumerate(seq_lines):
        if ensembl_protein and ensembl_protein in line:
            human_idx = i
            break
        if "homo_sapiens" in line.lower() or "human" in line.lower():
            human_idx = i
            break

    if human_idx is not None and human_idx > 0:
        human_line = seq_lines.pop(human_idx)
        seq_lines.insert(0, human_line)
        phylip_path.write_text(header + "\n" + "\n".join(seq_lines) + "\n")


# =====================================================================
#  ENTRY POINT
# =====================================================================


def main():
    """Run pervasive selection analysis from the command line."""
    ap = argparse.ArgumentParser(
        description="Pervasive positive selection analysis (CodeML site models)"
    )
    ap.add_argument("--ids-json", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument(
        "--msa", required=True, help="Protein MSA from the evolution module (orthologs.fasta)"
    )
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    with open(args.ids_json) as f:
        ids = json.load(f)
    record = SeqIO.read(args.fasta, "fasta")
    sequence = str(record.seq)
    msa_path = Path(args.msa)
    outdir = Path(args.outdir)

    # Handle missing MSA (sentinel file from Nextflow)
    if not msa_path.exists() or msa_path.stat().st_size == 0:
        log.warning("No MSA available (%s) – writing empty pervasive selection TSV.", msa_path)
        rows = [{"position": i + 1, "aa": sequence[i]} for i in range(len(sequence))]
        outpath = outdir / "pervasive_selection.tsv"
        with open(outpath, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        log.info("Wrote %s (empty – no MSA)", outpath)
        return

    rows = annotate_pervasive_selection(ids, sequence, msa_path, outdir)

    outpath = outdir / "pervasive_selection.tsv"
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    log.info("Wrote %s", outpath)


if __name__ == "__main__":
    main()
