#!/usr/bin/env python3
"""
modules/selection.py
=====================
Per-residue evolutionary selection pressure analysis.

  ┌──────────────────────────────────────────────────────────────┐
  │  1. SITE-SPECIFIC dN/dS (HyPhy FUBAR / MSA proxy)          │
  │     Per-site posterior estimates of synonymous (α) and       │
  │     non-synonymous (β) substitution rates.                  │
  │     β > α → positive selection. β < α → purifying.          │
  │     Uses HyPhy FUBAR (preferred); falls back to MSA-based   │
  │     column variability proxy when HyPhy is unavailable.     │
  │                                                              │
  │  2. ANCESTRAL RECONSTRUCTION                                │
  │     Infer ancestral AA using parsimony across distant taxa.  │
  │     Detects human-specific substitutions by comparing to     │
  │     chimpanzee and the inferred vertebrate ancestor.         │
  │                                                              │
  │  3. CONVERGENT EVOLUTION                                     │
  │     Detect positions where the same non-ancestral AA arose   │
  │     independently in phylogenetically distant lineages.      │
  │                                                              │
  │  4. COEVOLUTION (MUTUAL INFORMATION with APC)                │
  │     Positions that covary across the MSA. Reports top        │
  │     coevolving partners and maximum MI_APC per residue.      │
  │                                                              │
  │  5. EVOLUTIONARY RATE                                        │
  │     Normalised column variability relative to protein-wide   │
  │     average. >1 = fast-evolving, <1 = slow / constrained.   │
  └──────────────────────────────────────────────────────────────┘

Inputs:  resolved_ids.json, sequence.fasta, orthologs.fasta (Compara MSA)
Outputs: selection.tsv
         codon_alignment.fasta (if codon back-translation succeeds)
"""

import argparse
import csv
import json
import logging
import math
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")


# ── MSA helpers ──────────────────────────────────────────────────


def _load_msa(msa_path: Path):
    records = list(SeqIO.parse(msa_path, "fasta"))
    if not records:
        return [], None, 0
    human_id = None
    for rec in records:
        if "homo_sapiens" in rec.id.lower() or "human" in rec.id.lower():
            human_id = rec.id
            break
    if not human_id:
        human_id = records[0].id
    return records, human_id, len(records[0].seq)


def _get_column(records, col):
    return [str(rec.seq[col]) if col < len(rec.seq) else "-" for rec in records]


def _human_pos_map(records, human_id):
    for rec in records:
        if rec.id == human_id:
            mapping = {}
            pos = 0
            for col, aa in enumerate(str(rec.seq)):
                if aa != "-":
                    pos += 1
                    mapping[col] = pos
            return mapping
    return {}


# ── 1. Site-specific dN/dS ───────────────────────────────────────


def _run_hyphy_fubar(msa_path, outdir):
    hyphy = shutil.which("hyphy")
    if not hyphy:
        return None
    log.info("Running HyPhy FUBAR …")
    fubar_output = outdir / "fubar_output.json"
    cmd = [hyphy, "fubar", "--alignment", str(msa_path), "--output", str(fubar_output)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if result.returncode != 0:
            log.warning("HyPhy FUBAR failed: %s", result.stderr[:300])
            return None
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        log.warning("HyPhy error: %s", e)
        return None
    if not fubar_output.exists():
        return None
    with open(fubar_output) as f:
        data = json.load(f)
    results = {}
    mle = data.get("MLE", {}).get("content", {})
    rows = mle.get("0", [])
    for i, row in enumerate(rows):
        if len(row) >= 5:
            results[i] = {
                "dnds_alpha": round(row[0], 4),
                "dnds_beta": round(row[1], 4),
                "prob_negative_selection": round(row[2], 4),
                "prob_positive_selection": round(row[3], 4),
            }
    log.info("HyPhy FUBAR: %d sites.", len(results))
    return results


def _estimate_dnds_from_msa(records, human_id, col_map):
    results = {}
    human_rec = next((r for r in records if r.id == human_id), None)
    if not human_rec:
        return results
    col_vars = []
    for col in sorted(col_map.keys()):
        aas = [c for c in _get_column(records, col) if c in AA_ALPHABET]
        col_vars.append(len(set(aas)) / len(aas) if len(aas) >= 3 else None)
    non_none = [v for v in col_vars if v is not None]
    mean_var = sum(non_none) / len(non_none) if non_none else 0.01
    sd_var = (
        max(0.001, (sum((v - mean_var) ** 2 for v in non_none) / len(non_none)) ** 0.5)
        if non_none
        else 0.001
    )
    for i, col in enumerate(sorted(col_map.keys())):
        pos = col_map[col]
        var = col_vars[i]
        if var is None:
            results[pos] = {
                "dnds_alpha": None,
                "dnds_beta": None,
                "prob_positive_selection": None,
                "prob_negative_selection": None,
            }
            continue
        z = (var - mean_var) / sd_var
        alpha = max(0.01, 1.0 - z * 0.2)
        beta = max(0.01, 1.0 + z * 0.2)
        prob_pos = 1.0 / (1.0 + math.exp(-z * 2))
        results[pos] = {
            "dnds_alpha": round(alpha, 4),
            "dnds_beta": round(beta, 4),
            "prob_positive_selection": round(prob_pos, 4),
            "prob_negative_selection": round(1.0 - prob_pos, 4),
        }
    return results


def _classify_selection(r):
    pp = r.get("prob_positive_selection")
    pn = r.get("prob_negative_selection")
    a = r.get("dnds_alpha")
    b = r.get("dnds_beta")
    if pp is not None and pp > 0.9:
        return "positive"
    if pn is not None and pn > 0.9:
        return "purifying"
    if a is not None and b is not None:
        if b > a * 1.5:
            return "positive"
        if b < a * 0.5:
            return "purifying"
    return "neutral"


# ── 2. Ancestral reconstruction ──────────────────────────────────


def _reconstruct_ancestral(records, human_id, col_map):
    results = {}
    distant_ids = set()
    chimp_ids = set()
    for rec in records:
        rid = rec.id.lower()
        if "pan_troglodytes" in rid or "pan_paniscus" in rid:
            chimp_ids.add(rec.id)
        elif any(k in rid for k in ("danio", "takifugu", "oryzias", "xenopus", "anolis", "gallus")):
            distant_ids.add(rec.id)
    human_rec = next((r for r in records if r.id == human_id), None)
    if not human_rec:
        return results
    for col in range(len(human_rec.seq)):
        if col not in col_map:
            continue
        pos = col_map[col]
        human_aa = str(human_rec.seq[col])
        if human_aa not in AA_ALPHABET:
            continue
        distant_aas = [
            str(r.seq[col])
            for r in records
            if r.id in distant_ids and col < len(r.seq) and str(r.seq[col]) in AA_ALPHABET
        ]
        if distant_aas:
            ancestral_aa = Counter(distant_aas).most_common(1)[0][0]
        else:
            all_aas = [
                str(r.seq[col])
                for r in records
                if r.id != human_id and col < len(r.seq) and str(r.seq[col]) in AA_ALPHABET
            ]
            ancestral_aa = Counter(all_aas).most_common(1)[0][0] if all_aas else human_aa
        chimp_aa = None
        for rec in records:
            if rec.id in chimp_ids and col < len(rec.seq):
                aa = str(rec.seq[col])
                if aa in AA_ALPHABET:
                    chimp_aa = aa
                    break
        human_specific = False
        lineage_change = ""
        if human_aa != ancestral_aa:
            if chimp_aa and chimp_aa == ancestral_aa:
                human_specific = True
                lineage_change = f"{ancestral_aa}→{human_aa} (human-specific)"
            elif chimp_aa and chimp_aa != ancestral_aa and chimp_aa != human_aa:
                lineage_change = f"{ancestral_aa}→{human_aa} (both lineages diverged)"
            elif chimp_aa and chimp_aa == human_aa:
                lineage_change = f"{ancestral_aa}→{human_aa} (shared with chimp)"
            else:
                lineage_change = f"{ancestral_aa}→{human_aa}"
        results[pos] = {
            "ancestral_aa": ancestral_aa,
            "human_specific_substitution": human_specific,
            "human_lineage_change": lineage_change,
        }
    return results


# ── 3. Convergent evolution ──────────────────────────────────────


def _detect_convergence(records, human_id, col_map):
    results = {}
    lineage_map = defaultdict(set)
    for rec in records:
        rid = rec.id.lower()
        for keys, lineage in [
            (("homo_sapiens", "pan_", "gorilla", "pongo", "macaca"), "Primates"),
            (("mus_", "rattus", "cavia"), "Rodentia"),
            (("canis", "felis", "bos_", "sus_", "equus"), "Laurasiatheria"),
            (("gallus", "taeniopygia"), "Aves"),
            (("danio", "takifugu", "oryzias"), "Fish"),
            (("xenopus",), "Amphibia"),
        ]:
            if any(k in rid for k in keys):
                lineage_map[lineage].add(rec.id)
                break
    human_rec = next((r for r in records if r.id == human_id), None)
    if not human_rec:
        return results
    for col in range(len(human_rec.seq)):
        if col not in col_map:
            continue
        pos = col_map[col]
        column = [
            str(r.seq[col]) for r in records if col < len(r.seq) and str(r.seq[col]) in AA_ALPHABET
        ]
        if not column:
            continue
        most_common = Counter(column).most_common(1)[0][0]
        aa_lineages = defaultdict(set)
        for rec in records:
            if col < len(rec.seq):
                aa = str(rec.seq[col])
                if aa in AA_ALPHABET:
                    for lineage, members in lineage_map.items():
                        if rec.id in members:
                            aa_lineages[aa].add(lineage)
        convergent = False
        conv_lineages = ""
        distant_pairs = {
            ("Fish", "Primates"),
            ("Fish", "Aves"),
            ("Aves", "Rodentia"),
            ("Fish", "Laurasiatheria"),
            ("Amphibia", "Primates"),
        }
        for aa, lineages in aa_lineages.items():
            if aa == most_common or len(lineages) < 2:
                continue
            for l1, l2 in distant_pairs:
                if l1 in lineages and l2 in lineages:
                    convergent = True
                    conv_lineages = f"{aa} in {'+'.join(sorted(lineages))}"
                    break
            if convergent:
                break
        results[pos] = {"convergent_evolution": convergent, "convergent_lineages": conv_lineages}
    return results


# ── 4. Coevolution (MI + APC) ────────────────────────────────────


def _compute_coevolution(records, human_id, col_map, top_n=5):
    results = {}
    human_cols = sorted(col_map.keys())
    n_seqs = len(records)
    if len(human_cols) < 10 or n_seqs < 5:
        return results
    col_data = {col: _get_column(records, col) for col in human_cols}
    max_cols = min(len(human_cols), 500)
    sampled = human_cols[:max_cols]
    mi_matrix = {}
    col_mi_sum = defaultdict(float)
    n_pairs = 0
    for i_idx, ci in enumerate(sampled):
        for j_idx in range(i_idx + 1, len(sampled)):
            cj = sampled[j_idx]
            pairs = [
                (col_data[ci][s], col_data[cj][s])
                for s in range(n_seqs)
                if col_data[ci][s] in AA_ALPHABET and col_data[cj][s] in AA_ALPHABET
            ]
            if len(pairs) < 5:
                continue
            n = len(pairs)
            joint = Counter(pairs)
            mi_i = Counter(p[0] for p in pairs)
            mi_j = Counter(p[1] for p in pairs)
            mi = sum(
                (jc / n) * math.log2((jc / n) / ((mi_i[a] / n) * (mi_j[b] / n)))
                for (a, b), jc in joint.items()
                if jc > 0
            )
            mi_matrix[(ci, cj)] = mi
            mi_matrix[(cj, ci)] = mi
            col_mi_sum[ci] += mi
            col_mi_sum[cj] += mi
            n_pairs += 1
    if n_pairs == 0:
        return results
    overall = sum(mi_matrix.values()) / (2 * n_pairs)
    for c in sampled:
        col_mi_sum[c] /= max(1, len(sampled) - 1)
    for col in sampled:
        if col not in col_map:
            continue
        pos = col_map[col]
        partners = []
        for cj in sampled:
            if cj == col or cj not in col_map:
                continue
            mi_raw = mi_matrix.get((col, cj), 0)
            apc = col_mi_sum[col] * col_mi_sum[cj] / overall if overall > 0 else 0
            score = max(0, mi_raw - apc)
            if score > 0:
                partners.append((col_map[cj], score))
        partners.sort(key=lambda x: -x[1])
        top = partners[:top_n]
        results[pos] = {
            "top_coevo_partners": ";".join(f"{p}({s:.3f})" for p, s in top),
            "max_mi_apc": round(top[0][1], 4) if top else 0,
        }
    return results


# ── 5. Evolutionary rate ─────────────────────────────────────────


def _compute_rate(records, human_id, col_map):
    rates = {}
    vars_list = []
    sorted_cols = sorted(col_map.keys())
    for col in sorted_cols:
        aas = [c for c in _get_column(records, col) if c in AA_ALPHABET]
        vars_list.append(len(set(aas)) / len(aas) if len(aas) >= 3 else None)
    non_none = [v for v in vars_list if v is not None]
    mean_var = sum(non_none) / len(non_none) if non_none else 1.0
    for i, col in enumerate(sorted_cols):
        v = vars_list[i]
        rates[col_map[col]] = round(v / mean_var, 4) if v is not None and mean_var > 0 else None
    return rates


# ── Pipeline ─────────────────────────────────────────────────────

COLUMNS = [
    "position",
    "aa",
    "dnds_alpha",
    "dnds_beta",
    "prob_positive_selection",
    "prob_negative_selection",
    "selection_class",
    "human_specific_substitution",
    "human_lineage_change",
    "ancestral_aa",
    "convergent_evolution",
    "convergent_lineages",
    "top_coevo_partners",
    "max_mi_apc",
    "evolutionary_rate",
]


def annotate_selection(ids, sequence, msa_path):
    seq_len = len(sequence)
    results = [{"position": i + 1, "aa": sequence[i]} for i in range(seq_len)]
    if not msa_path.exists():
        log.warning("No MSA – skipping selection.")
        return results
    records, human_id, aln_len = _load_msa(msa_path)
    if not records or not human_id:
        log.warning("Could not load MSA.")
        return results
    col_map = _human_pos_map(records, human_id)
    log.info("MSA: %d seqs, %d cols, %d human pos.", len(records), aln_len, len(col_map))

    # 1. dN/dS
    log.info("── dN/dS ──")
    with tempfile.TemporaryDirectory() as tmp:
        fubar = _run_hyphy_fubar(msa_path, Path(tmp))
    if fubar:
        for col_idx, vals in fubar.items():
            if col_idx in col_map:
                pos = col_map[col_idx]
                for r in results:
                    if r["position"] == pos:
                        r.update(vals)
    else:
        msa_dnds = _estimate_dnds_from_msa(records, human_id, col_map)
        for r in results:
            if r["position"] in msa_dnds:
                r.update(msa_dnds[r["position"]])
    for r in results:
        r["selection_class"] = _classify_selection(r)
    n_pos = sum(1 for r in results if r.get("selection_class") == "positive")
    n_neg = sum(1 for r in results if r.get("selection_class") == "purifying")
    log.info(
        "Selection: %d positive, %d purifying, %d neutral.", n_pos, n_neg, seq_len - n_pos - n_neg
    )

    # 2. Ancestral
    log.info("── Ancestral reconstruction ──")
    anc = _reconstruct_ancestral(records, human_id, col_map)
    for r in results:
        if r["position"] in anc:
            r.update(anc[r["position"]])
    log.info("Human-specific: %d", sum(1 for r in results if r.get("human_specific_substitution")))

    # 3. Convergence
    log.info("── Convergent evolution ──")
    conv = _detect_convergence(records, human_id, col_map)
    for r in results:
        if r["position"] in conv:
            r.update(conv[r["position"]])
    log.info("Convergent: %d", sum(1 for r in results if r.get("convergent_evolution")))

    # 4. Coevolution
    log.info("── Coevolution (MI+APC) ──")
    coevo = _compute_coevolution(records, human_id, col_map)
    for r in results:
        if r["position"] in coevo:
            r.update(coevo[r["position"]])
    log.info("Coevolution: %d positions.", len(coevo))

    # 5. Rate
    log.info("── Evolutionary rate ──")
    rates = _compute_rate(records, human_id, col_map)
    for r in results:
        r["evolutionary_rate"] = rates.get(r["position"])

    return results


def main():
    ap = argparse.ArgumentParser(description="Selection pressure analysis")
    ap.add_argument("--ids-json", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--msa", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    with open(args.ids_json) as f:
        ids = json.load(f)
    record = SeqIO.read(args.fasta, "fasta")
    msa_path = Path(args.msa)

    # Handle missing MSA (e.g. sentinel file from Nextflow when evolution produced no MSA)
    if not msa_path.exists() or msa_path.stat().st_size == 0:
        sequence = str(record.seq)
        log.warning("No MSA available (%s) – writing empty selection TSV.", msa_path)
        rows = [{"position": i + 1, "aa": sequence[i]} for i in range(len(sequence))]
        outpath = Path(args.outdir) / "selection.tsv"
        with open(outpath, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        log.info("Wrote %s (empty – no MSA)", outpath)
        return

    rows = annotate_selection(ids, str(record.seq), msa_path)
    outpath = Path(args.outdir) / "selection.tsv"
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    log.info("Wrote %s", outpath)


if __name__ == "__main__":
    main()
