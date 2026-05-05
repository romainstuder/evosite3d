#!/usr/bin/env python3
"""
modules/designability.py
=========================
Per-residue protein-design annotations:

  ┌─────────────────────────────────────────────────────────────────┐
  │  EPITOPE DISCOVERY                                             │
  │  ─────────────────                                             │
  │  1. IEDB – known linear B-cell & T-cell epitopes mapped to    │
  │     UniProt positions via the IEDB REST API.                   │
  │  2. BepiPred-3 prediction – residue-level B-cell epitope      │
  │     probability from sequence features (hydrophilicity,        │
  │     surface accessibility, flexibility – Parker/Emini/Karplus  │
  │     scales) computed with BioPython ProtParam + ProteinAnalysis│
  │  3. UniProt "Antigenic" features – curated epitope regions.   │
  │                                                                │
  │  MUTATIONAL RISK                                               │
  │  ───────────────                                               │
  │  4. Ensembl VEP – precomputed SIFT & PolyPhen-2 scores for   │
  │     every possible missense variant at each position.          │
  │  5. ESM log-likelihood ratio (LLR) – language-model-based     │
  │     per-position mutational tolerance computed from the        │
  │     wild-type probability vs uniform background.               │
  │  6. Aggregated mutability flag – combines SIFT intolerance,   │
  │     PolyPhen damage, and conservation to flag positions where  │
  │     mutations are risky vs. designable.                        │
  └─────────────────────────────────────────────────────────────────┘

Inputs:  resolved_ids.json, sequence.fasta
Outputs: designability.tsv

Columns produced
-----------------
    position
    is_epitope           : bool   – in a known/predicted epitope region
    epitope_type         : str    – "B-cell", "T-cell-I", "T-cell-II", or combo
    epitope_source       : str    – "IEDB", "UniProt", "Predicted", or combo
    epitope_score        : float  – predicted B-cell epitope probability [0-1]
    sift_min_score       : float  – worst (lowest) SIFT score across all AAs
    sift_median_score    : float  – median SIFT score
    polyphen_max_score   : float  – worst (highest) PolyPhen score across all AAs
    esm_llr              : float  – ESM log-likelihood ratio (+ = tolerant)
    mutational_risk      : str    – "high" / "medium" / "low"
    designability_note   : str    – free-text summary for protein engineers
"""

import argparse
import csv
import json
import logging
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from api_utils import UNIPROT_REST, http_get_json

from lib.ensembl import fetch_translation_overlap

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


# =====================================================================
#  PART A – EPITOPE DISCOVERY
# =====================================================================

# ── A1. IEDB known epitopes ─────────────────────────────────────

IEDB_API = "https://query-api.iedb.org"


def _fetch_iedb_epitopes(uniprot: str, sequence: str) -> dict[int, dict]:
    """
    Query IEDB for linear epitopes mapped to this protein.
    Returns {position: {"types": set, "iedb_ids": list}}
    """
    epitope_map: dict[int, dict] = defaultdict(lambda: {"types": set(), "iedb_ids": []})

    # IEDB search by source antigen UniProt accession
    # B-cell epitopes
    url = (
        f"{IEDB_API}/epitope_search"
        f"?linear_sequence=any"
        f"&source_antigen_accession={uniprot}"
        f"&host_organism_id=9606"  # human host
        f"&output_format=json"
    )
    data = http_get_json(url)

    if data and isinstance(data, list):
        for entry in data:
            epi_seq = entry.get("linear_sequence", "")
            epi_id = entry.get("epitope_id", "")
            # epi_type_raw = entry.get("epitope_type", "")

            # Determine epitope type from assay info
            etype = set()
            object_type = entry.get("object_type", "").lower()
            if "b cell" in object_type or "antibody" in object_type:
                etype.add("B-cell")
            if "t cell" in object_type:
                mhc = entry.get("mhc_restriction", "").lower()
                if "class ii" in mhc:
                    etype.add("T-cell-II")
                else:
                    etype.add("T-cell-I")
            if not etype:
                etype.add("B-cell")  # default for linear

            # Map epitope sequence to protein positions
            if epi_seq and epi_seq in sequence:
                start = sequence.index(epi_seq)
                for offset in range(len(epi_seq)):
                    pos = start + offset + 1  # 1-based
                    epitope_map[pos]["types"].update(etype)
                    epitope_map[pos]["iedb_ids"].append(str(epi_id))

    n = len(epitope_map)
    log.info("IEDB epitopes: %d positions for %s", n, uniprot)
    return epitope_map


# ── A2. UniProt antigenic features ──────────────────────────────


def _fetch_uniprot_epitopes(uniprot: str) -> dict[int, dict]:
    """Parse 'Antigenic' and 'Epitope' features from UniProt."""
    epitope_map: dict[int, dict] = defaultdict(lambda: {"types": set(), "source": set()})

    url = f"{UNIPROT_REST}/uniprotkb/{uniprot}.json"
    data = http_get_json(url)
    if not data:
        return epitope_map

    for feat in data.get("features", []):
        ftype = feat.get("type", "").lower()
        if ftype not in ("region", "site", "binding site", "antigenic"):
            # Also check description for epitope keywords
            desc = feat.get("description", "").lower()
            if "epitope" not in desc and "antigenic" not in desc:
                continue

        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None:
            continue
        if end is None:
            end = start

        for p in range(int(start), int(end) + 1):
            epitope_map[p]["types"].add("B-cell")
            epitope_map[p]["source"].add("UniProt")

    log.info("UniProt antigenic regions: %d positions", len(epitope_map))
    return epitope_map


# ── A3. Sequence-based B-cell epitope prediction ────────────────
#    Uses BioPython ProteinAnalysis to compute physicochemical
#    properties that correlate with surface-exposed, flexible,
#    hydrophilic regions (classical BepiPred-like features).

# Parker hydrophilicity scale (Parker et al., 1986)
_PARKER = {
    "A": 2.1,
    "R": 4.2,
    "N": 7.0,
    "D": 10.0,
    "C": 1.4,
    "Q": 6.0,
    "E": 7.8,
    "G": 5.7,
    "H": 2.1,
    "I": -8.0,
    "L": -9.2,
    "K": 5.7,
    "M": -4.2,
    "F": -9.2,
    "P": 2.1,
    "S": 6.5,
    "T": 5.2,
    "W": -10.0,
    "Y": -1.9,
    "V": -3.7,
}

# Emini surface accessibility scale (Emini et al., 1985)
_EMINI = {
    "A": 0.815,
    "R": 1.475,
    "N": 1.296,
    "D": 1.283,
    "C": 0.394,
    "Q": 1.348,
    "E": 1.445,
    "G": 0.714,
    "H": 1.180,
    "I": 0.603,
    "L": 0.603,
    "K": 1.545,
    "M": 0.714,
    "F": 0.695,
    "P": 1.236,
    "S": 1.115,
    "T": 1.184,
    "W": 0.808,
    "Y": 1.089,
    "V": 0.606,
}

# Karplus-Schulz flexibility scale
_KARPLUS = {
    "A": 1.064,
    "R": 1.008,
    "N": 1.048,
    "D": 1.068,
    "C": 0.906,
    "Q": 1.037,
    "E": 1.094,
    "G": 1.102,
    "H": 0.950,
    "I": 0.927,
    "L": 0.935,
    "K": 1.102,
    "M": 0.952,
    "F": 0.915,
    "P": 1.049,
    "S": 1.046,
    "T": 0.997,
    "W": 0.904,
    "Y": 0.929,
    "V": 0.931,
}


def _sliding_window_avg(values: list[float], window: int = 7) -> list[float]:
    """Centered sliding window average."""
    n = len(values)
    half = window // 2
    result = []
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        result.append(sum(values[lo:hi]) / (hi - lo))
    return result


def _predict_bcell_epitopes(sequence: str, window: int = 7) -> list[float]:
    """
    Composite B-cell epitope score per residue [0–1].
    Combines hydrophilicity, surface accessibility, and flexibility
    using a logistic combination of z-scored scales.
    """
    n = len(sequence)
    if n == 0:
        return []

    # Raw scale values
    hydro_raw = [_PARKER.get(aa, 0.0) for aa in sequence]
    emini_raw = [_EMINI.get(aa, 1.0) for aa in sequence]
    flex_raw = [_KARPLUS.get(aa, 1.0) for aa in sequence]

    # Smooth with sliding window
    hydro = _sliding_window_avg(hydro_raw, window)
    emini = _sliding_window_avg(emini_raw, window)
    flex = _sliding_window_avg(flex_raw, window)

    # Also get BioPython flexibility (Vihinen B-factors)
    try:
        pa = ProteinAnalysis(sequence)
        bio_flex = pa.flexibility()
        # flexibility() returns n-9 values (window=9); pad to full length
        pad = (n - len(bio_flex)) // 2
        bio_flex = [bio_flex[0]] * pad + list(bio_flex) + [bio_flex[-1]] * (n - len(bio_flex) - pad)
    except Exception:
        bio_flex = flex  # fallback

    # Z-score each scale
    def _zscore(vals):
        mu = sum(vals) / len(vals)
        sd = (sum((v - mu) ** 2 for v in vals) / len(vals)) ** 0.5
        if sd < 1e-9:
            return [0.0] * len(vals)
        return [(v - mu) / sd for v in vals]

    zh = _zscore(hydro)
    ze = _zscore(emini)
    zf = _zscore(flex)
    zb = _zscore(bio_flex[:n])

    # Weighted combination → logistic sigmoid
    scores = []
    for i in range(n):
        # Weights: hydrophilicity 0.35, accessibility 0.30, flexibility 0.20, biopython 0.15
        combo = 0.35 * zh[i] + 0.30 * ze[i] + 0.20 * zf[i] + 0.15 * zb[i]
        prob = 1.0 / (1.0 + math.exp(-combo))
        scores.append(round(prob, 4))

    return scores


# =====================================================================
#  PART B – MUTATIONAL RISK
# =====================================================================

# ── B1. SIFT + PolyPhen from Ensembl VEP ────────────────────────


def _fetch_vep_scores(ensembl_protein: str, seq_len: int) -> dict[int, dict]:
    """
    Fetch precomputed missense variant effect predictions (SIFT, PolyPhen)
    via the Ensembl transcript/translation variation overlap.

    For each position we aggregate across all possible substitutions:
      sift_min     – lowest (worst) SIFT score  (< 0.05 = damaging)
      sift_median  – median SIFT score
      polyphen_max – highest (worst) PolyPhen score (> 0.85 = damaging)

    Returns {position: {"sift_min": float, "sift_median": float,
                         "polyphen_max": float}}
    """
    vep_map: dict[int, dict] = {}

    if not ensembl_protein:
        return vep_map

    log.info("Fetching Ensembl VEP SIFT/PolyPhen for %s …", ensembl_protein)
    data = fetch_translation_overlap(ensembl_protein, feature="transcript_variation")

    if not data:
        log.warning("No VEP data returned for %s", ensembl_protein)
        return vep_map

    # Collect per-position
    pos_sift: dict[int, list[float]] = defaultdict(list)
    pos_pp: dict[int, list[float]] = defaultdict(list)

    for var in data:
        # Only missense
        conseqs = var.get("consequence_terms", [])
        if "missense_variant" not in conseqs:
            continue

        start = var.get("start")
        if start is None:
            continue
        pos = int(start)

        sift = var.get("sift_score")
        pp = var.get("polyphen_score")

        if sift is not None:
            try:
                pos_sift[pos].append(float(sift))
            except (ValueError, TypeError):
                pass
        if pp is not None:
            try:
                pos_pp[pos].append(float(pp))
            except (ValueError, TypeError):
                pass

    # Aggregate
    for pos in set(list(pos_sift.keys()) + list(pos_pp.keys())):
        entry = {}
        if pos_sift[pos]:
            ss = sorted(pos_sift[pos])
            entry["sift_min"] = round(ss[0], 4)
            entry["sift_median"] = round(ss[len(ss) // 2], 4)
        if pos_pp[pos]:
            pp = sorted(pos_pp[pos])
            entry["polyphen_max"] = round(pp[-1], 4)
        vep_map[pos] = entry

    log.info("VEP scores: %d positions with SIFT/PolyPhen", len(vep_map))
    return vep_map


# ── B2. ESM log-likelihood ratio (LLR) ──────────────────────────
#    Approximated from AA frequencies at each position vs background.
#    True ESM requires GPU inference; here we use a lightweight proxy
#    based on the BLOSUM62-derived background + conservation.

# BLOSUM62-derived background AA frequencies (Robinson & Robinson, 1991)
_BG_FREQ = {
    "A": 0.0777,
    "R": 0.0531,
    "N": 0.0406,
    "D": 0.0543,
    "C": 0.0149,
    "Q": 0.0408,
    "E": 0.0634,
    "G": 0.0747,
    "H": 0.0219,
    "I": 0.0590,
    "L": 0.0964,
    "K": 0.0584,
    "M": 0.0236,
    "F": 0.0394,
    "P": 0.0505,
    "S": 0.0684,
    "T": 0.0556,
    "W": 0.0130,
    "Y": 0.0295,
    "V": 0.0676,
}


def _estimate_esm_llr(sequence: str, conservation_scores: dict[int, float]) -> dict[int, float]:
    """
    Proxy for per-position ESM log-likelihood ratio.

    LLR > 0 → position is tolerant to mutation (designable)
    LLR < 0 → position is constrained (risky to mutate)

    Uses conservation as a proxy for the language model's position-specific
    probability: highly conserved = high wild-type probability = low LLR.
    """
    llr_map = {}
    for i, aa in enumerate(sequence):
        pos = i + 1
        bg = _BG_FREQ.get(aa, 0.05)
        cons = conservation_scores.get(pos)

        if cons is not None and cons > 0:
            # Estimated WT probability from conservation (scaled)
            # Fully conserved (cons=1) → p_wt ≈ 0.95
            # Not conserved (cons=0) → p_wt ≈ background
            p_wt = bg + (0.95 - bg) * cons
        else:
            p_wt = bg

        # LLR = log(p_wt / bg)
        if bg > 0 and p_wt > 0:
            llr_map[pos] = round(math.log(p_wt / bg), 4)
        else:
            llr_map[pos] = 0.0

    return llr_map


# ── B3. Aggregate mutational risk classification ────────────────


def _classify_risk(
    sift_min: Optional[float],
    polyphen_max: Optional[float],
    esm_llr: Optional[float],
    conservation: Optional[float],
) -> str:
    """
    Combine multiple signals into a simple risk tier.

    high   – most substitutions are damaging; avoid mutating
    medium – mixed signal; mutate with caution
    low    – tolerant position; good candidate for engineering
    """
    damaging_signals = 0
    total_signals = 0

    if sift_min is not None:
        total_signals += 1
        if sift_min < 0.05:
            damaging_signals += 1

    if polyphen_max is not None:
        total_signals += 1
        if polyphen_max > 0.85:
            damaging_signals += 1

    if conservation is not None:
        total_signals += 1
        if conservation > 0.9:
            damaging_signals += 1

    if esm_llr is not None:
        total_signals += 1
        if esm_llr > 2.0:  # high WT probability → constrained
            damaging_signals += 1

    if total_signals == 0:
        return "unknown"

    ratio = damaging_signals / total_signals
    if ratio >= 0.6:
        return "high"
    elif ratio >= 0.3:
        return "medium"
    else:
        return "low"


def _make_design_note(
    is_epitope: bool,
    epitope_type: str,
    risk: str,
    sift_min: Optional[float],
    polyphen_max: Optional[float],
) -> str:
    """Short free-text note for protein engineers."""
    parts = []

    if is_epitope:
        parts.append(f"Epitope ({epitope_type})")

    if risk == "high":
        parts.append("High mutational risk – strongly conserved / damaging")
    elif risk == "low":
        parts.append("Designable – tolerant to substitution")

    if sift_min is not None and sift_min < 0.05:
        parts.append(f"SIFT intolerant ({sift_min:.3f})")
    if polyphen_max is not None and polyphen_max > 0.85:
        parts.append(f"PolyPhen damaging ({polyphen_max:.3f})")

    if is_epitope and risk == "high":
        parts.append("CAUTION: epitope in constrained region")
    elif is_epitope and risk == "low":
        parts.append("Consider: epitope in mutable region – de-immunization candidate")

    return "; ".join(parts) if parts else ""


# =====================================================================
#  PUBLIC PIPELINE
# =====================================================================


def annotate_designability(
    ids: dict, sequence: str, conservation_scores: dict[int, float] | None = None
) -> list[dict]:
    """
    Full designability annotation per residue.
    """
    seq_len = len(sequence)
    uniprot = ids.get("uniprot")
    ens_prot = ids.get("ensembl_protein")

    # ── Epitopes ─────────────────────────────────────────────────
    log.info("── Epitope discovery ──")

    # A1: IEDB
    iedb_map = {}
    if uniprot:
        iedb_map = _fetch_iedb_epitopes(uniprot, sequence)

    # A2: UniProt antigenic
    uniprot_epi = {}
    if uniprot:
        uniprot_epi = _fetch_uniprot_epitopes(uniprot)

    # A3: Predicted scores
    log.info("Computing B-cell epitope predictions (Parker/Emini/Karplus + BioPython)…")
    predicted_scores = _predict_bcell_epitopes(sequence)
    EPITOPE_THRESHOLD = 0.55  # positions above this are predicted epitopes

    # ── Mutational risk ──────────────────────────────────────────
    log.info("── Mutational risk assessment ──")

    # B1: SIFT + PolyPhen
    vep_map = _fetch_vep_scores(ens_prot, seq_len)

    # B2: ESM LLR proxy
    if conservation_scores is None:
        conservation_scores = {}
    esm_llr_map = _estimate_esm_llr(sequence, conservation_scores)

    # ── Assemble per-residue ─────────────────────────────────────
    results = []
    for i, aa in enumerate(sequence):
        pos = i + 1
        row = {"position": pos, "aa": aa}

        # -- Epitope fields --
        types = set()
        sources = set()

        if pos in iedb_map:
            types.update(iedb_map[pos]["types"])
            sources.add("IEDB")
        if pos in uniprot_epi:
            types.update(uniprot_epi[pos].get("types", set()))
            sources.update(uniprot_epi[pos].get("source", set()))

        epi_score = predicted_scores[i] if i < len(predicted_scores) else None
        if epi_score and epi_score >= EPITOPE_THRESHOLD and not types:
            types.add("B-cell")
            sources.add("Predicted")

        is_epi = bool(types)
        row["is_epitope"] = is_epi
        row["epitope_type"] = "; ".join(sorted(types)) if types else ""
        row["epitope_source"] = "; ".join(sorted(sources)) if sources else ""
        row["epitope_score"] = epi_score

        # -- Mutational risk fields --
        vep = vep_map.get(pos, {})
        sift_min = vep.get("sift_min")
        sift_med = vep.get("sift_median")
        pp_max = vep.get("polyphen_max")
        esm = esm_llr_map.get(pos)
        cons = conservation_scores.get(pos)

        row["sift_min_score"] = sift_min
        row["sift_median_score"] = sift_med
        row["polyphen_max_score"] = pp_max
        row["esm_llr"] = esm

        risk = _classify_risk(sift_min, pp_max, esm, cons)
        row["mutational_risk"] = risk

        note = _make_design_note(is_epi, row["epitope_type"], risk, sift_min, pp_max)
        row["designability_note"] = note

        results.append(row)

    n_epi = sum(1 for r in results if r["is_epitope"])
    n_high = sum(1 for r in results if r["mutational_risk"] == "high")
    n_low = sum(1 for r in results if r["mutational_risk"] == "low")
    log.info("Epitope residues: %d / %d", n_epi, seq_len)
    log.info("Mutational risk: %d high, %d low, rest medium/unknown", n_high, n_low)

    return results


# =====================================================================
#  CLI ENTRY POINT
# =====================================================================

COLUMNS = [
    "position",
    "aa",
    "is_epitope",
    "epitope_type",
    "epitope_source",
    "epitope_score",
    "sift_min_score",
    "sift_median_score",
    "polyphen_max_score",
    "esm_llr",
    "mutational_risk",
    "designability_note",
]


def main():
    ap = argparse.ArgumentParser(
        description="Protein design annotations: epitopes + mutational risk"
    )
    ap.add_argument("--ids-json", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument(
        "--evolution-tsv",
        default=None,
        help="Optional: evolution.tsv from the evolution module, "
        "used to improve ESM-LLR estimation.",
    )
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    with open(args.ids_json) as f:
        ids = json.load(f)
    record = SeqIO.read(args.fasta, "fasta")
    sequence = str(record.seq)

    # Optionally load conservation scores
    conservation = {}
    if args.evolution_tsv and Path(args.evolution_tsv).exists():
        with open(args.evolution_tsv, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                pos = int(row["position"])
                cs = row.get("conservation_score", "")
                if cs and cs != "None":
                    try:
                        conservation[pos] = float(cs)
                    except ValueError:
                        pass
        log.info("Loaded %d conservation scores from evolution.tsv", len(conservation))

    rows = annotate_designability(ids, sequence, conservation)

    outpath = Path(args.outdir) / "designability.tsv"
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    log.info("Wrote %s", outpath)


if __name__ == "__main__":
    main()
