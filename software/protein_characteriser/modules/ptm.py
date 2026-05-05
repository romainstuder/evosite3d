#!/usr/bin/env python3
"""
modules/ptm.py
===============
Per-residue post-translational modification annotations.

Sources:
  1. UniProt – curated "Modified residue", "Glycosylation",
     "Lipidation", "Cross-link", "Disulfide bond" features.
  2. dbPTM   – supplementary PTM database (best-effort).

Inputs:  resolved_ids.json, sequence.fasta
Outputs: ptm.tsv  (position, ptm)
"""

import argparse
import csv
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from api_utils import UNIPROT_REST, http_get_json

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

_PTM_FEATURE_TYPES = {
    "Modified residue",
    "Glycosylation",
    "Lipidation",
    "Cross-link",
    "Disulfide bond",
}


def fetch_uniprot_ptms(uniprot: str) -> dict[int, list[str]]:
    ptms: dict[int, list[str]] = defaultdict(list)
    url = f"{UNIPROT_REST}/uniprotkb/{uniprot}.json"
    data = http_get_json(url)
    if not data:
        return ptms

    for feat in data.get("features", []):
        ftype = feat.get("type", "")
        if ftype not in _PTM_FEATURE_TYPES:
            continue

        desc = (feat.get("description") or ftype).strip()
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None:
            continue

        if ftype == "Disulfide bond":
            ptms[int(start)].append("Disulfide bond")
            if end and end != start:
                ptms[int(end)].append("Disulfide bond")
        elif start == end or end is None:
            ptms[int(start)].append(desc)
        else:
            for p in range(int(start), int(end) + 1):
                ptms[p].append(desc)

    log.info("UniProt PTMs: %d positions for %s", len(ptms), uniprot)
    return ptms


def fetch_dbptm(uniprot: str) -> dict[int, list[str]]:
    ptms: dict[int, list[str]] = defaultdict(list)
    url = f"https://awi.cuhk.edu.cn/dbPTM/rest/ptm_annotation/{uniprot}"
    data = http_get_json(url)
    if not data or not isinstance(data, list):
        return ptms
    for entry in data:
        pos = entry.get("position")
        ptype = entry.get("ptm_type", "Unknown PTM")
        if pos:
            try:
                ptms[int(pos)].append(ptype)
            except (ValueError, TypeError):
                pass
    if ptms:
        log.info("dbPTM: %d additional positions for %s", len(ptms), uniprot)
    return ptms


def annotate_ptm(uniprot: str, sequence: str) -> list[dict]:
    seq_len = len(sequence)
    results = [{"position": i + 1, "aa": sequence[i], "ptm": ""} for i in range(seq_len)]
    if not uniprot:
        log.warning("No UniProt – skipping PTM.")
        return results

    ptm_map = fetch_uniprot_ptms(uniprot)

    # Supplementary (best-effort)
    try:
        dbptm = fetch_dbptm(uniprot)
        for pos, types in dbptm.items():
            existing = set(ptm_map.get(pos, []))
            for t in types:
                if t not in existing:
                    ptm_map[pos].append(t)
    except Exception as e:
        log.debug("dbPTM failed (non-critical): %s", e)

    for r in results:
        pos = r["position"]
        if pos in ptm_map:
            # Deduplicate preserving order
            seen = set()
            unique = []
            for p in ptm_map[pos]:
                if p not in seen:
                    seen.add(p)
                    unique.append(p)
            r["ptm"] = "; ".join(unique)

    n = sum(1 for r in results if r["ptm"])
    log.info("PTM: %d / %d residues annotated.", n, seq_len)
    return results


# ── Entry point ──────────────────────────────────────────────────


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids-json", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    with open(args.ids_json) as f:
        ids = json.load(f)
    record = SeqIO.read(args.fasta, "fasta")

    rows = annotate_ptm(ids.get("uniprot"), str(record.seq))

    outpath = Path(args.outdir) / "ptm.tsv"
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["position", "aa", "ptm"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    log.info("Wrote %s", outpath)


if __name__ == "__main__":
    main()
