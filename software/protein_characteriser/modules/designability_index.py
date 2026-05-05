#!/usr/bin/env python3
"""
modules/designability_index.py
===============================
Composite designability scoring function that integrates all annotation
signals into a single per-residue Designability Index (DI).

  ┌──────────────────────────────────────────────────────────────┐
  │  DESIGNABILITY INDEX (DI)                                   │
  │  ────────────────────────                                    │
  │  A weighted composite score [0–1] where:                    │
  │    1.0 = ideal for engineering (mutable, tolerant, exposed)  │
  │    0.0 = do not mutate (conserved, buried, functional)      │
  │                                                              │
  │  SIGNAL WEIGHTS (benchmarked on DMS datasets):              │
  │                                                              │
  │  Constraint signals (push DI toward 0):                     │
  │    - ESM-LLR (0.25)        : language model constraint      │
  │    - Conservation (0.15)    : evolutionary constraint        │
  │    - SIFT (0.10)           : sequence-based intolerance      │
  │    - PolyPhen (0.10)       : structure-based damage          │
  │    - Functional site (0.10) : active/binding/catalytic       │
  │                                                              │
  │  Opportunity signals (push DI toward 1):                    │
  │    - RSA/accessibility (0.10) : surface exposure             │
  │    - Disorder (0.08)        : flexibility / tolerance        │
  │    - Entropy (0.07)         : ESM position entropy           │
  │    - Non-interface (0.05)   : not at PPI interface           │
  │                                                              │
  │  Penalty flags (hard modifiers):                            │
  │    - Disease variant at position → DI × 0.3                 │
  │    - Cancer hotspot → DI × 0.2                              │
  │    - Active site / metal binding → DI × 0.1                 │
  │    - Disulfide bond → DI × 0.1                              │
  │    - Epitope → DI × 0.7 (mild penalty for immunogenicity)   │
  │                                                              │
  │  OUTPUT:                                                     │
  │    designability_index : float [0–1]                        │
  │    di_class : "engineer" / "caution" / "conserve" / "avoid" │
  │    di_rationale : human-readable explanation                 │
  │    di_signals : JSON of individual signal contributions     │
  └──────────────────────────────────────────────────────────────┘

This module reads the outputs of ALL other modules (via the merged TSV
or individual module TSVs) and computes the composite index.

Inputs:  All module TSVs (evolution, structure, annotation, biophysics,
         ptm, disease, cancer, designability, esm_scoring, af_multimer)
         OR the merged *_residues.tsv
Outputs: designability_index.tsv
"""

import argparse
import csv
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


# =====================================================================
#  Signal extraction functions
# =====================================================================


def _safe_float(val, default=None):
    if val is None or val == "" or val == "None":
        return default
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def _safe_bool(val):
    if isinstance(val, bool):
        return val
    if isinstance(val, str):
        return val.lower() in ("true", "1", "yes")
    return bool(val) if val else False


def _normalise(val, lo, hi, invert=False):
    """Normalise value to [0, 1]. If invert, 1 means low original value."""
    if val is None:
        return None
    clamped = max(lo, min(hi, val))
    norm = (clamped - lo) / (hi - lo) if hi > lo else 0.5
    return (1.0 - norm) if invert else norm


# =====================================================================
#  Composite scoring
# =====================================================================

# Weights for each signal (sum to 1.0)
WEIGHTS = {
    # Constraint signals (high value → constrained → low DI)
    "esm_llr": 0.25,
    "conservation": 0.15,
    "sift": 0.10,
    "polyphen": 0.10,
    "functional_site": 0.10,
    # Opportunity signals (high value → designable → high DI)
    "accessibility": 0.10,
    "disorder": 0.08,
    "esm_entropy": 0.07,
    "non_interface": 0.05,
}


def compute_designability_index(row: dict) -> dict:
    """
    Compute the composite Designability Index for a single residue.

    Parameters
    ----------
    row : dict with all annotation fields

    Returns
    -------
    dict with: designability_index, di_class, di_rationale, di_signals
    """
    signals = {}
    rationale_parts = []
    n_available = 0
    weighted_sum = 0.0
    weight_sum = 0.0

    # ── ESM LLR (constraint) ────────────────────────────────
    esm_llr = _safe_float(row.get("esm_llr"))
    if esm_llr is not None:
        # LLR > 0 means constrained; normalise: high LLR → high constraint → low DI
        # Typical range: -2 to +5
        s = _normalise(esm_llr, -2.0, 5.0, invert=True)
        signals["esm_llr"] = round(s, 4)
        weighted_sum += s * WEIGHTS["esm_llr"]
        weight_sum += WEIGHTS["esm_llr"]
        n_available += 1
        if esm_llr > 3.0:
            rationale_parts.append("strongly constrained by ESM")
        elif esm_llr < 0:
            rationale_parts.append("ESM-tolerant")

    # ── Conservation (constraint) ────────────────────────────
    cons = _safe_float(row.get("conservation_score"))
    if cons is not None:
        s = 1.0 - cons  # high conservation → low designability
        signals["conservation"] = round(s, 4)
        weighted_sum += s * WEIGHTS["conservation"]
        weight_sum += WEIGHTS["conservation"]
        n_available += 1
        if cons > 0.95:
            rationale_parts.append("near-invariant across vertebrates")

    # ── SIFT (constraint) ───────────────────────────────────
    sift = _safe_float(row.get("sift_min_score"))
    if sift is not None:
        # SIFT < 0.05 = intolerant; high SIFT = tolerant
        s = sift  # already 0–1 where 1 = tolerant
        signals["sift"] = round(s, 4)
        weighted_sum += s * WEIGHTS["sift"]
        weight_sum += WEIGHTS["sift"]
        n_available += 1

    # ── PolyPhen (constraint) ───────────────────────────────
    pp = _safe_float(row.get("polyphen_max_score"))
    if pp is not None:
        # PolyPhen > 0.85 = damaging; invert
        s = 1.0 - pp
        signals["polyphen"] = round(s, 4)
        weighted_sum += s * WEIGHTS["polyphen"]
        weight_sum += WEIGHTS["polyphen"]
        n_available += 1

    # ── Functional site (constraint) ─────────────────────────
    func = row.get("functional_site", "")
    has_func = bool(func and str(func).strip())
    s = 0.0 if has_func else 1.0
    signals["functional_site"] = s
    weighted_sum += s * WEIGHTS["functional_site"]
    weight_sum += WEIGHTS["functional_site"]
    n_available += 1
    if has_func:
        rationale_parts.append(f"functional site ({str(func)[:40]})")

    # ── Accessibility (opportunity) ──────────────────────────
    rsa = _safe_float(row.get("rsa"))
    if rsa is not None:
        s = _normalise(rsa, 0.0, 0.6)  # RSA > 0.6 considered fully exposed
        signals["accessibility"] = round(s, 4)
        weighted_sum += s * WEIGHTS["accessibility"]
        weight_sum += WEIGHTS["accessibility"]
        n_available += 1
    else:
        acc = row.get("accessibility", "")
        if acc == "exposed":
            s = 1.0
        elif acc == "intermediate":
            s = 0.5
        elif acc == "buried":
            s = 0.1
        else:
            s = None
        if s is not None:
            signals["accessibility"] = s
            weighted_sum += s * WEIGHTS["accessibility"]
            weight_sum += WEIGHTS["accessibility"]
            n_available += 1

    # ── Disorder (opportunity) ───────────────────────────────
    disorder = _safe_float(row.get("disorder_probability"))
    if disorder is not None:
        signals["disorder"] = round(disorder, 4)
        weighted_sum += disorder * WEIGHTS["disorder"]
        weight_sum += WEIGHTS["disorder"]
        n_available += 1

    # ── ESM entropy (opportunity) ────────────────────────────
    entropy_norm = _safe_float(row.get("esm_entropy_norm"))
    if entropy_norm is not None:
        signals["esm_entropy"] = round(entropy_norm, 4)
        weighted_sum += entropy_norm * WEIGHTS["esm_entropy"]
        weight_sum += WEIGHTS["esm_entropy"]
        n_available += 1

    # ── Non-interface (opportunity) ──────────────────────────
    in_interface = _safe_bool(row.get("in_interface"))
    predicted_interface = _safe_bool(row.get("predicted_interface"))
    is_interface = in_interface or predicted_interface
    s = 0.0 if is_interface else 1.0
    signals["non_interface"] = s
    weighted_sum += s * WEIGHTS["non_interface"]
    weight_sum += WEIGHTS["non_interface"]
    n_available += 1
    if is_interface:
        rationale_parts.append("at protein-protein interface")

    # ── Compute raw DI ───────────────────────────────────────
    if weight_sum > 0:
        raw_di = weighted_sum / weight_sum
    else:
        raw_di = 0.5  # no data

    # ── Hard penalty flags ───────────────────────────────────
    penalties = []

    if row.get("disease_association") and str(row["disease_association"]).strip():
        raw_di *= 0.3
        penalties.append("disease variant (×0.3)")

    if (
        _safe_bool(row.get("is_cancer_hotspot"))
        or (_safe_float(row.get("cosmic_mutation_count")) or 0) >= 5
    ):
        raw_di *= 0.2
        penalties.append("cancer hotspot (×0.2)")

    ptm = row.get("ptm", "")
    if ptm and "disulfide" in str(ptm).lower():
        raw_di *= 0.1
        penalties.append("disulfide bond (×0.1)")

    func_str = str(func).lower() if func else ""
    if "active site" in func_str or "catalytic" in func_str or "metal binding" in func_str:
        raw_di *= 0.1
        penalties.append("catalytic/metal site (×0.1)")

    if _safe_bool(row.get("is_epitope")):
        raw_di *= 0.7
        penalties.append("epitope (×0.7)")

    final_di = round(max(0.0, min(1.0, raw_di)), 4)

    # ── Classification ───────────────────────────────────────
    if final_di >= 0.7:
        di_class = "engineer"
    elif final_di >= 0.4:
        di_class = "caution"
    elif final_di >= 0.15:
        di_class = "conserve"
    else:
        di_class = "avoid"

    # ── Rationale ────────────────────────────────────────────
    if penalties:
        rationale_parts.extend(penalties)

    if not rationale_parts:
        if di_class == "engineer":
            rationale_parts.append("no constraint signals – good candidate")
        elif di_class == "avoid":
            rationale_parts.append("multiple constraint signals")

    return {
        "designability_index": final_di,
        "di_class": di_class,
        "di_rationale": "; ".join(rationale_parts),
        "di_signals": json.dumps(signals),
        "di_n_signals": n_available,
    }


# =====================================================================
#  ENTRY POINT
# =====================================================================

COLUMNS = [
    "position",
    "aa",
    "designability_index",
    "di_class",
    "di_rationale",
    "di_signals",
    "di_n_signals",
]


def main():
    ap = argparse.ArgumentParser(
        description="Compute composite Designability Index from merged annotations"
    )
    ap.add_argument(
        "--merged-tsv",
        required=True,
        help="Path to the merged *_residues.tsv (with all annotations)",
    )
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    # Load merged TSV
    rows = []
    with open(args.merged_tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)

    log.info("Computing Designability Index for %d residues …", len(rows))

    results = []
    for row in rows:
        pos = row.get("position", "")
        di = compute_designability_index(row)
        di["position"] = pos
        di["aa"] = row.get("amino_acid", "")
        results.append(di)

    # Stats
    classes = defaultdict(int)
    dis = []
    for r in results:
        classes[r["di_class"]] += 1
        dis.append(r["designability_index"])

    avg_di = sum(dis) / len(dis) if dis else 0
    log.info(
        "DI distribution: engineer=%d  caution=%d  conserve=%d  avoid=%d  (avg=%.3f)",
        classes["engineer"],
        classes["caution"],
        classes["conserve"],
        classes["avoid"],
        avg_di,
    )

    outpath = Path(args.outdir) / "designability_index.tsv"
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(results)
    log.info("Wrote %s", outpath)


if __name__ == "__main__":
    main()
