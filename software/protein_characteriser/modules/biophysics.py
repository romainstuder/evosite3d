#!/usr/bin/env python3
"""
modules/biophysics.py
======================
Per-residue biophysical, functional site, and sequence-context annotations.

  ┌──────────────────────────────────────────────────────────────┐
  │  1. PHYSICOCHEMICAL PROPERTIES (BioPython ProtParam)         │
  │     Per-residue charge, hydrophobicity (Kyte-Doolittle),     │
  │     molecular weight, polarity class, amino acid class.      │
  │                                                              │
  │  2. INTRINSIC DISORDER (MobiDB-lite / IUPred proxy)          │
  │     Disorder probability from MobiDB REST API, supplemented  │
  │     by a sequence-composition proxy when unavailable.        │
  │                                                              │
  │  3. FUNCTIONAL SITES (UniProt)                               │
  │     Active sites, metal-binding, catalytic residues,         │
  │     nucleotide/DNA/RNA-binding residues, calcium-binding.    │
  │                                                              │
  │  4. SHORT LINEAR MOTIFS (ELM)                                │
  │     Eukaryotic Linear Motif instances from the ELM DB,       │
  │     providing SLiM-level functional annotation.              │
  │                                                              │
  │  5. ISOFORM COVERAGE                                         │
  │     Whether the residue is present in all known isoforms     │
  │     or is alternatively spliced (from UniProt isoforms).     │
  │                                                              │
  │  6. CODON CONTEXT                                            │
  │     CpG context and rare codon flag from the Ensembl CDS,   │
  │     relevant for expression optimisation and mutagenesis.    │
  │                                                              │
  │  7. PROTEOMICS EVIDENCE                                      │
  │     Peptide-level proteomics evidence from PeptideAtlas /    │
  │     neXtProt via UniProt, indicating experimental detection. │
  └──────────────────────────────────────────────────────────────┘

Inputs:  resolved_ids.json, sequence.fasta
Outputs: biophysics.tsv
"""

import argparse
import csv
import json
import logging
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from api_utils import UNIPROT_REST, http_get_json

from lib.ensembl import fetch_cds_for_protein

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


# =====================================================================
#  1. PHYSICOCHEMICAL PROPERTIES
# =====================================================================

# Kyte-Doolittle hydrophobicity scale
_KD_HYDROPHOBICITY = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}

# Charge at pH 7.0
_CHARGE_PH7 = {
    "R": +1,
    "K": +1,
    "H": +0.1,  # His is ~10% protonated at pH 7
    "D": -1,
    "E": -1,
}

# Amino acid classification
_AA_CLASS = {
    "G": "tiny",
    "A": "tiny",
    "V": "hydrophobic",
    "I": "hydrophobic",
    "L": "hydrophobic",
    "M": "hydrophobic",
    "F": "aromatic",
    "W": "aromatic",
    "Y": "aromatic",
    "P": "special",
    "C": "special",
    "S": "polar",
    "T": "polar",
    "N": "polar",
    "Q": "polar",
    "H": "charged",
    "K": "charged",
    "R": "charged",
    "D": "charged",
    "E": "charged",
}

# Residue molecular weights (monoisotopic, Da)
_MW = {
    "A": 71.04,
    "R": 156.10,
    "N": 114.04,
    "D": 115.03,
    "C": 103.01,
    "Q": 128.06,
    "E": 129.04,
    "G": 57.02,
    "H": 137.06,
    "I": 113.08,
    "L": 113.08,
    "K": 128.09,
    "M": 131.04,
    "F": 147.07,
    "P": 97.05,
    "S": 87.03,
    "T": 101.05,
    "W": 186.08,
    "Y": 163.06,
    "V": 99.07,
}

# Volume (Å³, Zamyatin 1972)
_VOLUME = {
    "A": 88.6,
    "R": 173.4,
    "N": 114.1,
    "D": 111.1,
    "C": 108.5,
    "Q": 143.8,
    "E": 138.4,
    "G": 60.1,
    "H": 153.2,
    "I": 166.7,
    "L": 166.7,
    "K": 168.6,
    "M": 162.9,
    "F": 189.9,
    "P": 112.7,
    "S": 89.0,
    "T": 116.1,
    "W": 227.8,
    "Y": 193.6,
    "V": 140.0,
}


def compute_physicochemical(sequence: str) -> list[dict]:
    """Per-residue physicochemical properties."""
    results = []

    # Windowed hydrophobicity (window=9)
    window = 9
    half = window // 2
    n = len(sequence)

    for i, aa in enumerate(sequence):
        row = {}
        row["hydrophobicity"] = _KD_HYDROPHOBICITY.get(aa, 0.0)
        row["charge_ph7"] = _CHARGE_PH7.get(aa, 0.0)
        row["residue_mw"] = _MW.get(aa, 0.0)
        row["residue_volume"] = _VOLUME.get(aa, 0.0)
        row["aa_class"] = _AA_CLASS.get(aa, "unknown")

        # Windowed hydrophobicity (GRAVY-like)
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        window_aas = sequence[lo:hi]
        row["hydrophobicity_window9"] = round(
            sum(_KD_HYDROPHOBICITY.get(a, 0) for a in window_aas) / len(window_aas), 3
        )

        results.append(row)

    return results


# =====================================================================
#  2. INTRINSIC DISORDER (MobiDB)
# =====================================================================


def fetch_disorder(uniprot: str, seq_len: int) -> dict[int, dict]:
    """
    Fetch per-residue disorder annotations from MobiDB.
    Returns {position: {"disorder_probability": float,
                         "is_disordered": bool,
                         "disorder_source": str}}
    """
    disorder: dict[int, dict] = {}

    url = f"https://mobidb.bio.unipd.it/api/download?acc={uniprot}&format=json"
    data = http_get_json(url)

    if data and isinstance(data, dict):
        # MobiDB returns consensus disorder regions
        consensus = data.get("consensus", {})
        disorder_regions = consensus.get("disorder", {}).get("regions", [])

        # Mark disordered positions
        disordered_set = set()
        for region in disorder_regions:
            if isinstance(region, list) and len(region) >= 2:
                for p in range(int(region[0]), int(region[1]) + 1):
                    disordered_set.add(p)

        # Also get prediction scores if available
        prediction = data.get("prediction-disorder-mobidb_lite", {})
        scores = prediction.get("scores", [])

        for pos in range(1, seq_len + 1):
            score = None
            if scores and pos <= len(scores):
                try:
                    score = float(scores[pos - 1])
                except (ValueError, TypeError, IndexError):
                    pass

            is_dis = pos in disordered_set
            disorder[pos] = {
                "disorder_probability": round(score, 4) if score is not None else "",
                "is_disordered": is_dis,
                "disorder_source": "MobiDB" if (score is not None or is_dis) else "",
            }

        log.info("MobiDB disorder: %d / %d positions disordered", len(disordered_set), seq_len)

    # Fallback: composition-based disorder proxy using BioPython
    if not disorder:
        log.info("MobiDB unavailable – using composition-based disorder proxy.")
        disorder = _predict_disorder_proxy(uniprot, seq_len)

    return disorder


def _predict_disorder_proxy(uniprot: str, seq_len: int) -> dict[int, dict]:
    """
    Lightweight disorder estimation based on amino acid composition
    in a sliding window. Disorder-promoting residues: A, R, G, Q, S, P, E, K.
    Order-promoting: W, C, F, I, Y, V, L, N.
    """
    # Get sequence
    url = f"{UNIPROT_REST}/uniprotkb/{uniprot}.json"
    data = http_get_json(url)
    if not data:
        return {}
    sequence = data.get("sequence", {}).get("value", "")
    if not sequence:
        return {}

    _DISORDER_PROPENSITY = {
        "A": 0.06,
        "R": 0.18,
        "N": -0.01,
        "D": 0.19,
        "C": -0.20,
        "Q": 0.16,
        "E": 0.24,
        "G": 0.17,
        "H": -0.03,
        "I": -0.49,
        "L": -0.34,
        "K": 0.20,
        "M": -0.19,
        "F": -0.42,
        "P": 0.41,
        "S": 0.14,
        "T": -0.05,
        "W": -0.49,
        "Y": -0.34,
        "V": -0.39,
    }

    window = 21
    half = window // 2
    n = len(sequence)
    disorder = {}

    for i in range(n):
        pos = i + 1
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        win_seq = sequence[lo:hi]
        score = sum(_DISORDER_PROPENSITY.get(aa, 0) for aa in win_seq) / len(win_seq)
        # Sigmoid transform
        prob = 1.0 / (1.0 + math.exp(-10 * score))
        disorder[pos] = {
            "disorder_probability": round(prob, 4),
            "is_disordered": prob > 0.5,
            "disorder_source": "Composition-proxy",
        }

    return disorder


# =====================================================================
#  3. FUNCTIONAL SITES (UniProt features)
# =====================================================================


def fetch_functional_sites(uniprot: str) -> dict[int, list[str]]:
    """
    Active sites, metal-binding, catalytic residues, nucleotide/DNA/RNA
    binding, calcium binding from UniProt features.
    """
    sites: dict[int, list[str]] = defaultdict(list)

    url = f"{UNIPROT_REST}/uniprotkb/{uniprot}.json"
    data = http_get_json(url)
    if not data:
        return sites

    _SITE_TYPES = {
        "Active site",
        "Binding site",
        "Metal binding",
        "Site",
        "Calcium binding",
        "Nucleotide binding",
        "DNA binding",
    }

    for feat in data.get("features", []):
        ftype = feat.get("type", "")
        if ftype not in _SITE_TYPES:
            continue

        desc = feat.get("description", "")
        label = f"{ftype}"
        if desc:
            label += f": {desc}"

        # Get ligand info if available
        ligand = feat.get("ligand", {})
        if ligand:
            lig_name = ligand.get("name", "")
            if lig_name:
                label += f" [{lig_name}]"

        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None:
            continue
        if end is None:
            end = start

        for p in range(int(start), int(end) + 1):
            sites[p].append(label)

    log.info("Functional sites: %d positions for %s", len(sites), uniprot)
    return sites


# =====================================================================
#  4. SHORT LINEAR MOTIFS (ELM)
# =====================================================================


def fetch_elm_motifs(uniprot: str, sequence: str) -> dict[int, list[str]]:
    """
    Fetch ELM (Eukaryotic Linear Motif) instances for this protein.
    Falls back to regex-based scanning of known ELM patterns.
    """
    motifs: dict[int, list[str]] = defaultdict(list)

    # Try ELM REST API
    url = f"http://elm.eu.org/api/search/{uniprot}.json"
    data = http_get_json(url)

    if data and isinstance(data, dict):
        instances = data.get("instances", [])
        for inst in instances:
            elm_id = inst.get("elm_identifier", "")
            elm_name = inst.get("elm_name", elm_id)
            start = inst.get("start")
            end = inst.get("end")
            if start and end:
                label = f"{elm_id}: {elm_name}" if elm_name != elm_id else elm_id
                for p in range(int(start), int(end) + 1):
                    motifs[p].append(label)

    # If ELM API is unavailable, scan for a few high-value motifs
    if not motifs:
        log.info("ELM API unavailable – scanning for common motifs.")
        _COMMON_MOTIFS = {
            "NxS/T glycosylation": r"N[^P][ST]",
            "RGD cell attachment": r"RGD",
            "PEST degradation": r"[PE]{2,}[ST]",
            "KEN box (APC/C)": r"KEN",
            "D-box (APC/C)": r"R..L",
            "Nuclear export (NES)": r"L...L..LI",
            "SUMO conjugation": r"[VILMAFP]K.E",
            "SH3 binding": r"P..P",
            "WW domain binding": r"PP.Y",
            "14-3-3 binding": r"R[SFYW].[ST]",
        }
        for motif_name, pattern in _COMMON_MOTIFS.items():
            for match in re.finditer(pattern, sequence):
                start = match.start() + 1  # 1-based
                end = match.end()
                for p in range(start, end + 1):
                    motifs[p].append(motif_name)

    log.info("ELM/motifs: %d positions annotated", len(motifs))
    return motifs


# =====================================================================
#  5. ISOFORM COVERAGE
# =====================================================================


def fetch_isoform_coverage(uniprot: str, seq_len: int) -> dict[int, dict]:
    """
    Check if each residue is present in all UniProt isoforms
    or is alternatively spliced.
    """
    coverage: dict[int, dict] = {}

    url = f"{UNIPROT_REST}/uniprotkb/{uniprot}.json"
    data = http_get_json(url)
    if not data:
        return coverage

    # Count isoforms
    comments = data.get("comments", [])
    n_isoforms = 1
    for comment in comments:
        if comment.get("commentType") == "ALTERNATIVE PRODUCTS":
            isoforms = comment.get("isoforms", [])
            n_isoforms = max(n_isoforms, len(isoforms))

    # Find splice variants from features
    splice_positions = set()
    for feat in data.get("features", []):
        ftype = feat.get("type", "")
        if ftype not in ("Alternative sequence", "Splice variant"):
            continue
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start:
            if end is None:
                end = start
            for p in range(int(start), int(end) + 1):
                splice_positions.add(p)

    for pos in range(1, seq_len + 1):
        is_spliced = pos in splice_positions
        coverage[pos] = {
            "n_isoforms": n_isoforms,
            "is_alternatively_spliced": is_spliced,
            "isoform_coverage": "partial" if is_spliced else "all",
        }

    log.info(
        "Isoforms: %d total, %d alternatively spliced positions", n_isoforms, len(splice_positions)
    )
    return coverage


# =====================================================================
#  6. CODON CONTEXT
# =====================================================================

# Human codon usage (frequency per thousand, from Kazusa)
_RARE_CODONS = {
    # Codons used <10 per thousand in human
    "TCG",
    "CCG",
    "ACG",
    "GCG",  # xCG codons
    "CGA",
    "CGU",  # Arg rare codons
    "CUA",  # Leu rare
    "AUA",  # Ile rare
    "GUA",  # Val rare
}


def fetch_codon_context(ensembl_protein: str, sequence: str) -> dict[int, dict]:
    """
    Fetch CDS from Ensembl, extract per-residue codon, CpG status,
    and rare codon flag.
    """
    codons: dict[int, dict] = {}

    if not ensembl_protein:
        return codons

    cds = fetch_cds_for_protein(ensembl_protein)
    if not cds:
        return codons

    # Walk codons
    n_codons = len(cds) // 3
    for i in range(min(n_codons, len(sequence))):
        pos = i + 1
        codon = cds[i * 3 : i * 3 + 3].upper()
        if len(codon) < 3:
            break

        is_cpg = "CG" in codon
        is_rare = codon in _RARE_CODONS

        codons[pos] = {
            "codon": codon,
            "is_cpg_codon": is_cpg,
            "is_rare_codon": is_rare,
        }

    log.info(
        "Codon context: %d positions, %d CpG, %d rare",
        len(codons),
        sum(1 for c in codons.values() if c["is_cpg_codon"]),
        sum(1 for c in codons.values() if c["is_rare_codon"]),
    )
    return codons


# =====================================================================
#  7. PROTEOMICS EVIDENCE
# =====================================================================


def fetch_proteomics_evidence(uniprot: str) -> dict[int, str]:
    """
    Check neXtProt / UniProt for proteomics-level evidence.
    Returns {position: evidence_level} where evidence_level is
    "protein_level", "transcript_level", "homology", or "predicted".
    """
    evidence: dict[int, str] = {}

    # Global protein existence level from UniProt
    url = f"{UNIPROT_REST}/uniprotkb/{uniprot}.json"
    data = http_get_json(url)
    if not data:
        return evidence

    pe_level = data.get("proteinExistence", "")
    # Map to human-readable
    pe_map = {
        "1: Evidence at protein level": "protein_level",
        "2: Evidence at transcript level": "transcript_level",
        "3: Inferred from homology": "homology",
        "4: Predicted": "predicted",
        "5: Uncertain": "uncertain",
    }
    pe_str = pe_map.get(pe_level, pe_level)

    # Check for peptide-level features (mass spec evidence)
    peptide_positions = set()
    for feat in data.get("features", []):
        if feat.get("type") in ("Peptide", "Chain"):
            # These indicate experimentally observed peptides
            if any("mass spectrometry" in str(ev).lower() for ev in feat.get("evidences", [])):
                loc = feat.get("location", {})
                start = loc.get("start", {}).get("value")
                end = loc.get("end", {}).get("value")
                if start and end:
                    for p in range(int(start), int(end) + 1):
                        peptide_positions.add(p)

    seq_len = len(data.get("sequence", {}).get("value", ""))
    for pos in range(1, seq_len + 1):
        if pos in peptide_positions:
            evidence[pos] = "ms_detected"
        else:
            evidence[pos] = pe_str

    return evidence


# =====================================================================
#  PUBLIC PIPELINE
# =====================================================================

COLUMNS = [
    "position",
    "aa",
    # Physicochemical
    "hydrophobicity",
    "hydrophobicity_window9",
    "charge_ph7",
    "residue_mw",
    "residue_volume",
    "aa_class",
    # Disorder
    "disorder_probability",
    "is_disordered",
    # Functional sites
    "functional_site",
    # Linear motifs
    "elm_motif",
    # Isoform
    "n_isoforms",
    "is_alternatively_spliced",
    # Codon
    "codon",
    "is_cpg_codon",
    "is_rare_codon",
    # Proteomics
    "proteomics_evidence",
]


def annotate_biophysics(ids: dict, sequence: str) -> list[dict]:
    """Full biophysical / functional context annotation."""
    seq_len = len(sequence)
    uniprot = ids.get("uniprot")
    ens_prot = ids.get("ensembl_protein")

    # 1. Physicochemical
    log.info("Computing physicochemical properties…")
    physchem = compute_physicochemical(sequence)

    # 2. Disorder
    disorder = {}
    if uniprot:
        disorder = fetch_disorder(uniprot, seq_len)

    # 3. Functional sites
    sites = {}
    if uniprot:
        sites = fetch_functional_sites(uniprot)

    # 4. ELM motifs
    motifs = {}
    if uniprot:
        motifs = fetch_elm_motifs(uniprot, sequence)

    # 5. Isoform coverage
    isoforms = {}
    if uniprot:
        isoforms = fetch_isoform_coverage(uniprot, seq_len)

    # 6. Codon context
    codons = fetch_codon_context(ens_prot, sequence)

    # 7. Proteomics
    proteomics = {}
    if uniprot:
        proteomics = fetch_proteomics_evidence(uniprot)

    # ── Assemble ─────────────────────────────────────────────
    results = []
    for i in range(seq_len):
        pos = i + 1
        row = {"position": pos, "aa": sequence[i]}

        # Physicochemical
        pc = physchem[i] if i < len(physchem) else {}
        row["hydrophobicity"] = pc.get("hydrophobicity", "")
        row["hydrophobicity_window9"] = pc.get("hydrophobicity_window9", "")
        row["charge_ph7"] = pc.get("charge_ph7", 0.0)
        row["residue_mw"] = pc.get("residue_mw", "")
        row["residue_volume"] = pc.get("residue_volume", "")
        row["aa_class"] = pc.get("aa_class", "")

        # Disorder
        dis = disorder.get(pos, {})
        row["disorder_probability"] = dis.get("disorder_probability", "")
        row["is_disordered"] = dis.get("is_disordered", "")

        # Functional sites
        s = sites.get(pos, [])
        row["functional_site"] = "; ".join(s) if s else ""

        # ELM motifs
        m = motifs.get(pos, [])
        row["elm_motif"] = "; ".join(sorted(set(m))) if m else ""

        # Isoform
        iso = isoforms.get(pos, {})
        row["n_isoforms"] = iso.get("n_isoforms", "")
        row["is_alternatively_spliced"] = iso.get("is_alternatively_spliced", "")

        # Codon
        cod = codons.get(pos, {})
        row["codon"] = cod.get("codon", "")
        row["is_cpg_codon"] = cod.get("is_cpg_codon", "")
        row["is_rare_codon"] = cod.get("is_rare_codon", "")

        # Proteomics
        row["proteomics_evidence"] = proteomics.get(pos, "")

        results.append(row)

    log.info("Biophysics: %d residues annotated.", seq_len)
    return results


# ── Entry point ──────────────────────────────────────────────────


def main():
    ap = argparse.ArgumentParser(
        description="Biophysical, functional site, and sequence-context annotations"
    )
    ap.add_argument("--ids-json", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    with open(args.ids_json) as f:
        ids = json.load(f)
    record = SeqIO.read(args.fasta, "fasta")
    sequence = str(record.seq)

    rows = annotate_biophysics(ids, sequence)

    outpath = Path(args.outdir) / "biophysics.tsv"
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    log.info("Wrote %s", outpath)


if __name__ == "__main__":
    main()
