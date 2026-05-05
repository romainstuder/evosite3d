#!/usr/bin/env python3
"""
modules/esm_scoring.py
=======================
Per-residue protein language model scoring using ESM-2.

  ┌──────────────────────────────────────────────────────────────┐
  │  TRUE ESM-2 INFERENCE (replaces conservation-based proxy)   │
  │                                                              │
  │  1. Log-likelihood ratio (LLR) per position:                │
  │     LLR = log P(wt_aa | context) - log P_background(wt_aa)  │
  │     High LLR → position is strongly preferred by the model  │
  │     → constrained, risky to mutate                          │
  │     Low/negative LLR → position is flexible → designable    │
  │                                                              │
  │  2. Per-position entropy:                                    │
  │     Shannon entropy of the predicted AA distribution.        │
  │     Low entropy → constrained. High → tolerant.             │
  │                                                              │
  │  3. Full mutational landscape matrix (20×L):                │
  │     ΔLL for every possible substitution at every position.   │
  │     Exported as a separate matrix TSV for deep analysis.     │
  │                                                              │
  │  4. Top tolerated substitutions per position:                │
  │     The 3 best non-WT amino acids by ESM probability.       │
  │                                                              │
  │  Requires: pip install fair-esm torch                        │
  │  GPU recommended but CPU works for proteins < 1000 aa.      │
  │  Falls back to a BLOSUM62-based proxy if ESM unavailable.   │
  └──────────────────────────────────────────────────────────────┘

Inputs:  resolved_ids.json, sequence.fasta
Outputs: esm_scoring.tsv                 (per-residue summary)
         esm_mutation_matrix.tsv         (20×L substitution landscape)
         esm_jalview_features.txt        (Jalview features: constraint tiers)
         esm_jalview_annotations.txt     (Jalview bar graphs: LLR + entropy)
"""

import argparse
import csv
import logging
import math
import sys
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"

# BLOSUM62-derived background frequencies
_BG_FREQ = {
    "A": 0.0777,
    "C": 0.0149,
    "D": 0.0543,
    "E": 0.0634,
    "F": 0.0394,
    "G": 0.0747,
    "H": 0.0219,
    "I": 0.0590,
    "K": 0.0584,
    "L": 0.0964,
    "M": 0.0236,
    "N": 0.0406,
    "P": 0.0505,
    "Q": 0.0408,
    "R": 0.0531,
    "S": 0.0684,
    "T": 0.0556,
    "V": 0.0676,
    "W": 0.0130,
    "Y": 0.0295,
}


# =====================================================================
#  ESM-2 inference
# =====================================================================


def _try_esm_inference(sequence: str) -> dict | None:
    """
    Run ESM-2 (esm2_t33_650M_UR50D) inference on the sequence.

    Returns dict with:
        logits: list of dicts, one per position, {AA: log_prob}
        model_name: str
    Or None if ESM is not available.
    """
    try:
        import esm
        import torch
    except ImportError:
        return None

    log.info("Loading ESM-2 model (esm2_t33_650M_UR50D) …")
    try:
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    except Exception:
        # Try smaller model as fallback
        try:
            log.info("Falling back to esm2_t12_35M_UR50D …")
            model, alphabet = esm.pretrained.esm2_t12_35M_UR50D()
        except Exception as e:
            log.warning("ESM model loading failed: %s", e)
            return None

    batch_converter = alphabet.get_batch_converter()
    model.eval()

    # Use GPU if available
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    log.info("Running ESM-2 on %s (%d aa) …", device, len(sequence))

    data = [("protein", sequence)]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_tokens = batch_tokens.to(device)

    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[], return_contacts=False)

    # Extract log-probabilities per position
    # logits shape: (1, seq_len+2, vocab_size) — includes BOS/EOS tokens
    logits = results["logits"][0]  # (seq_len+2, vocab_size)

    # Convert to log-probabilities
    log_probs = torch.nn.functional.log_softmax(logits, dim=-1)

    # Map ESM alphabet tokens to standard AAs
    per_position = []
    for i in range(1, len(sequence) + 1):  # skip BOS token at index 0
        pos_probs = {}
        for aa in AA_ORDER:
            token_idx = alphabet.get_idx(aa)
            pos_probs[aa] = float(log_probs[i, token_idx].cpu())
        per_position.append(pos_probs)

    return {
        "logits": per_position,
        "model_name": "esm2_t33_650M_UR50D",
    }


def _blosum_proxy(sequence: str) -> dict:
    """
    BLOSUM62-based proxy when ESM is not available.
    Uses position-independent AA frequencies as a simple baseline.
    """
    log.info("Using BLOSUM62-based proxy (ESM not available).")
    per_position = []
    for aa in sequence:
        pos_probs = {}
        for target_aa in AA_ORDER:
            if target_aa == aa:
                pos_probs[target_aa] = math.log(0.7)  # pseudo WT preference
            else:
                pos_probs[target_aa] = math.log(_BG_FREQ.get(target_aa, 0.05))
        per_position.append(pos_probs)

    return {
        "logits": per_position,
        "model_name": "blosum62_proxy",
    }


# =====================================================================
#  Scoring functions
# =====================================================================


def compute_esm_scores(sequence: str) -> tuple[list[dict], list[dict]]:
    """
    Compute per-residue ESM scores and full mutation matrix.

    Returns:
        (per_residue_scores, mutation_matrix_rows)
    """
    # Try real ESM, fall back to proxy
    esm_result = _try_esm_inference(sequence)
    if esm_result is None:
        esm_result = _blosum_proxy(sequence)

    model_name = esm_result["model_name"]
    logits = esm_result["logits"]

    per_residue = []
    matrix_rows = []

    for i, aa in enumerate(sequence):
        pos = i + 1
        probs = logits[i]

        # ── LLR: log P(wt) - log P_background(wt) ───────────
        wt_log_prob = probs.get(aa, -5.0)
        bg_log_prob = math.log(_BG_FREQ.get(aa, 0.05))
        llr = round(wt_log_prob - bg_log_prob, 4)

        # ── Entropy of predicted distribution ────────────────
        # Convert log-probs to probs, compute Shannon entropy
        max_lp = max(probs.values())
        raw_probs = {a: math.exp(lp - max_lp) for a, lp in probs.items()}
        total = sum(raw_probs.values())
        norm_probs = {a: p / total for a, p in raw_probs.items()}

        entropy = -sum(p * math.log2(p) for p in norm_probs.values() if p > 1e-10)
        entropy = round(entropy, 4)
        max_entropy = math.log2(20)  # 4.322 bits

        # ── Top tolerated substitutions ──────────────────────
        # Sort non-WT AAs by probability
        non_wt = {a: norm_probs[a] for a in AA_ORDER if a != aa}
        top3 = sorted(non_wt.items(), key=lambda x: -x[1])[:3]
        top_subs = ";".join(f"{a}({p:.3f})" for a, p in top3)

        # ── Delta log-likelihood for each substitution ───────
        delta_row = {"position": pos, "wt_aa": aa}
        for target_aa in AA_ORDER:
            dll = round(probs.get(target_aa, -10) - wt_log_prob, 4)
            delta_row[f"dll_{target_aa}"] = dll
        matrix_rows.append(delta_row)

        per_residue.append(
            {
                "position": pos,
                "aa": aa,
                "esm_llr": llr,
                "esm_entropy": entropy,
                "esm_entropy_norm": round(entropy / max_entropy, 4),
                "esm_wt_prob": round(norm_probs.get(aa, 0), 4),
                "esm_top_substitutions": top_subs,
                "esm_model": model_name,
            }
        )

    log.info("ESM scoring complete (%s): %d positions.", model_name, len(per_residue))
    return per_residue, matrix_rows


# =====================================================================
#  Jalview visualisation tracks
# =====================================================================


def write_jalview_features(results: list[dict], seq_id: str, outpath: Path):
    """Write a Jalview features file with ESM constraint tiers and top substitutions.

    Produces colour-coded features for three constraint tiers based on
    the log-likelihood ratio (constrained LLR > 2, moderate 0–2,
    flexible < 0) and a separate group listing the top tolerated
    substitutions per position.

    Args:
        results: Per-residue score dicts as returned by
            :func:`compute_esm_scores`.
        seq_id: Sequence identifier to attach features to.
        outpath: Destination file path.
    """
    with open(outpath, "w") as f:
        # Colour definitions
        f.write("esm_constrained\t#d73027\n")
        f.write("esm_moderate\t#fee08b\n")
        f.write("esm_flexible\t#1a9850\n")

        # Constraint tier features
        f.write("\nSTARTGROUP\tESM Constraint\n")
        for r in results:
            pos = r["position"]
            llr = r["esm_llr"]
            if llr > 2.0:
                ftype = "esm_constrained"
                desc = f"Constrained (LLR={llr:.2f})"
            elif llr > 0.0:
                ftype = "esm_moderate"
                desc = f"Moderate (LLR={llr:.2f})"
            else:
                ftype = "esm_flexible"
                desc = f"Flexible (LLR={llr:.2f})"
            f.write(f"{desc}\t{seq_id}\t-1\t{pos}\t{pos}\t{ftype}\n")
        f.write("ENDGROUP\tESM Constraint\n")

        # Top substitutions
        f.write("\nSTARTGROUP\tESM Top Substitutions\n")
        for r in results:
            pos = r["position"]
            top_subs = r.get("esm_top_substitutions", "")
            if top_subs:
                f.write(f"Top subs: {top_subs}\t{seq_id}\t-1\t{pos}\t{pos}\tesm_substitution\n")
        f.write("ENDGROUP\tESM Top Substitutions\n")

    log.info("Wrote Jalview features: %s", outpath)


def write_jalview_annotations(results: list[dict], seq_id: str, outpath: Path):
    """Write a Jalview annotation file with LLR and entropy bar graphs.

    Produces two ``BAR_GRAPH`` tracks:

    * **ESM LLR** — log-likelihood ratio per position, coloured from
      green (flexible) through yellow to red (constrained).
    * **ESM Entropy** — normalised Shannon entropy, coloured from red
      (low entropy, constrained) to green (high entropy, tolerant).

    Args:
        results: Per-residue score dicts as returned by
            :func:`compute_esm_scores`.
        seq_id: Sequence identifier to attach annotations to.
        outpath: Destination file path.
    """
    with open(outpath, "w") as f:
        f.write("JALVIEW_ANNOTATION\n")
        f.write("# ESM-2 protein language model per-residue scores\n\n")
        f.write(f"SEQUENCE_REF\t{seq_id}\n")

        # LLR bar graph
        llr_vals = []
        for r in results:
            llr = r["esm_llr"]
            if llr > 2.0:
                colour = "d73027"
            elif llr > 1.0:
                colour = "ff8800"
            elif llr > 0.0:
                colour = "ffcc00"
            elif llr > -1.0:
                colour = "88cc00"
            else:
                colour = "1a9850"
            llr_vals.append(f"{llr:.3f},{llr:.2f},{colour}")

        f.write("BAR_GRAPH\tESM LLR\tLog-likelihood ratio vs background (high=constrained)\t")
        f.write("|".join(llr_vals))
        f.write("\n")

        # Entropy bar graph
        ent_vals = []
        for r in results:
            ent = r["esm_entropy_norm"]
            if ent > 0.7:
                colour = "1a9850"
            elif ent > 0.5:
                colour = "88cc00"
            elif ent > 0.3:
                colour = "ffcc00"
            elif ent > 0.15:
                colour = "ff8800"
            else:
                colour = "d73027"
            ent_vals.append(f"{ent:.3f},{ent:.2f},{colour}")

        f.write("BAR_GRAPH\tESM Entropy\tNormalised Shannon entropy (high=tolerant)\t")
        f.write("|".join(ent_vals))
        f.write("\n")

    log.info("Wrote Jalview annotations: %s", outpath)


# =====================================================================
#  ENTRY POINT
# =====================================================================

COLUMNS = [
    "position",
    "aa",
    "esm_llr",
    "esm_entropy",
    "esm_entropy_norm",
    "esm_wt_prob",
    "esm_top_substitutions",
    "esm_model",
]

MATRIX_COLUMNS = ["position", "wt_aa"] + [f"dll_{aa}" for aa in AA_ORDER]


def main():
    ap = argparse.ArgumentParser(description="ESM-2 protein language model scoring per residue")
    ap.add_argument("--ids-json", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    # with open(args.ids_json) as f:
    #     ids = json.load(f)
    record = SeqIO.read(args.fasta, "fasta")
    sequence = str(record.seq)
    seq_id = record.id
    outdir = Path(args.outdir)

    per_residue, matrix_rows = compute_esm_scores(sequence)

    # Per-residue summary
    with open(outdir / "esm_scoring.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(per_residue)
    log.info("Wrote esm_scoring.tsv")

    # Full mutation matrix (20×L)
    with open(outdir / "esm_mutation_matrix.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=MATRIX_COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(matrix_rows)
    log.info(
        "Wrote esm_mutation_matrix.tsv (%d positions × %d AAs)", len(matrix_rows), len(AA_ORDER)
    )

    # Jalview visualisation tracks
    write_jalview_features(per_residue, seq_id, outdir / "esm_jalview_features.txt")
    write_jalview_annotations(per_residue, seq_id, outdir / "esm_jalview_annotations.txt")


if __name__ == "__main__":
    main()
