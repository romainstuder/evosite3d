#!/usr/bin/env python3
"""
modules/cancer.py
==================
Per-residue somatic cancer mutation annotations.

  ┌──────────────────────────────────────────────────────────────┐
  │  SOMATIC MUTATIONS                                           │
  │  ─────────────────                                           │
  │  1. COSMIC (via Ensembl) – somatic variants overlapping the  │
  │     translation, filtered for pathogenic/COSMIC-annotated.   │
  │     Includes mutation count and COSMIC IDs.                  │
  │                                                              │
  │  2. cBioPortal cancer hotspots – recurrently mutated         │
  │     positions across pan-cancer datasets.                    │
  │                                                              │
  │  3. Cancer Gene Census context – whether the gene is in the  │
  │     COSMIC Cancer Gene Census (tumor suppressor / oncogene). │
  │                                                              │
  │  4. Aggregated cancer tier:                                  │
  │     • hotspot   – recurrently mutated (≥5 COSMIC samples)   │
  │     • moderate  – observed somatic mutations (2-4 samples)   │
  │     • rare      – ≤1 somatic observation                    │
  │     • none      – no somatic data                           │
  └──────────────────────────────────────────────────────────────┘

Data sources
------------
- Ensembl REST /overlap/translation – somatic variants with COSMIC xrefs
- cancerhotspots.org v2 API – statistically significant hotspots
- UniProt "Mutagenesis" features – curated functional mutagenesis data

Note: COSMIC's full database requires a license for direct download.
This module uses the publicly available COSMIC annotations exposed
through Ensembl and supplementary open databases.

Inputs:  resolved_ids.json, sequence.fasta
Outputs: cancer.tsv
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

from lib.ensembl import fetch_translation_overlap

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


# ── 1. COSMIC somatic variants via Ensembl ───────────────────────


def _fetch_cosmic_via_ensembl(ensembl_protein: str) -> dict[int, dict]:
    """
    Fetch somatic variants overlapping the translation from Ensembl.
    Ensembl integrates COSMIC data and exposes it via the overlap API.

    Returns {position: {"cosmic_ids": [...], "count": int,
                         "aa_changes": [...], "cancer_types": set()}}
    """
    somatic: dict[int, dict] = defaultdict(
        lambda: {"cosmic_ids": [], "count": 0, "aa_changes": [], "cancer_types": set()}
    )

    if not ensembl_protein:
        return somatic

    log.info("Fetching somatic variants from Ensembl for %s …", ensembl_protein)

    # Ensembl somatic variation overlay
    data = fetch_translation_overlap(ensembl_protein, feature="somatic_transcript_variation")

    if not data:
        # Fallback: try regular variation and filter for somatic
        data = fetch_translation_overlap(ensembl_protein, feature="transcript_variation")
        if data:
            # Filter for somatic / COSMIC
            data = [
                v
                for v in data
                if v.get("somatic")
                or any("COSM" in str(s) or "COSV" in str(s) for s in v.get("synonyms", []))
            ]
        else:
            data = []

    log.info("Raw somatic variants from Ensembl: %d", len(data) if data else 0)

    for var in data:
        start = var.get("start")
        if start is None:
            continue
        pos = int(start)

        var_id = var.get("id", "")
        residues = var.get("residues", "")
        conseqs = var.get("consequence_terms", [])

        # Only keep protein-affecting variants
        protein_affecting = {
            "missense_variant",
            "stop_gained",
            "frameshift_variant",
            "inframe_deletion",
            "inframe_insertion",
            "start_lost",
            "stop_lost",
            "splice_region_variant",
        }
        if not any(c in protein_affecting for c in conseqs):
            continue

        # Check for COSMIC ID
        is_cosmic = "COSM" in var_id or "COSV" in var_id
        synonyms = var.get("synonyms", [])
        cosmic_ids = [s for s in synonyms if "COSM" in str(s) or "COSV" in str(s)]
        if is_cosmic:
            cosmic_ids.append(var_id)

        # Deduplicate
        cosmic_ids = list(set(cosmic_ids))

        # Clinical significance can sometimes indicate cancer
        # clinical = var.get("clinical_significance", [])

        somatic[pos]["cosmic_ids"].extend(cosmic_ids)
        somatic[pos]["count"] += 1
        if residues:
            somatic[pos]["aa_changes"].append(residues)

        # Extract cancer type from phenotype if available
        phenotypes = var.get("phenotype", [])
        if isinstance(phenotypes, list):
            for pheno in phenotypes:
                if isinstance(pheno, dict):
                    desc = pheno.get("description", "")
                    if desc:
                        somatic[pos]["cancer_types"].add(desc)
                elif isinstance(pheno, str):
                    somatic[pos]["cancer_types"].add(pheno)

    # Deduplicate cosmic_ids per position
    for pos in somatic:
        somatic[pos]["cosmic_ids"] = list(set(somatic[pos]["cosmic_ids"]))

    log.info("COSMIC somatic positions: %d", len(somatic))
    return somatic


# ── 2. Cancer hotspots (cancerhotspots.org) ──────────────────────


def _fetch_cancer_hotspots(gene_symbol: str) -> dict[int, dict]:
    """
    Query cancerhotspots.org for statistically significant
    recurrently mutated residues.

    Returns {position: {"hotspot_score": float, "hotspot_type": str,
                         "cancer_types": set()}}
    """
    hotspots: dict[int, dict] = {}

    if not gene_symbol:
        return hotspots

    # cancerhotspots.org REST API
    url = f"https://www.cancerhotspots.org/api/hotspots/single/{gene_symbol}"
    data = http_get_json(url)

    if data and isinstance(data, list):
        for entry in data:
            pos = entry.get("residue")
            if pos is None:
                continue

            # Parse residue number from format like "R248"
            try:
                res_num = int("".join(c for c in str(pos) if c.isdigit()))
            except ValueError:
                continue

            hotspot_type = entry.get("type", "single")
            q_value = entry.get("qValue") or entry.get("q_value")
            score = None
            if q_value is not None:
                try:
                    score = float(q_value)
                except (ValueError, TypeError):
                    pass

            cancer_types = set()
            compositions = entry.get("cancerTypeComposition") or entry.get("composition", {})
            if isinstance(compositions, dict):
                cancer_types.update(compositions.keys())

            hotspots[res_num] = {
                "hotspot_score": score,
                "hotspot_type": hotspot_type,
                "cancer_types": cancer_types,
            }

    log.info("Cancer hotspots: %d positions for %s", len(hotspots), gene_symbol)
    return hotspots


# ── 3. UniProt mutagenesis data ──────────────────────────────────


def _fetch_uniprot_mutagenesis(uniprot: str) -> dict[int, list[str]]:
    """
    Parse UniProt "Mutagenesis" features – curated experimental
    mutagenesis results that often include cancer-relevant functional data.
    """
    muts: dict[int, list[str]] = defaultdict(list)

    url = f"{UNIPROT_REST}/uniprotkb/{uniprot}.json"
    data = http_get_json(url)
    if not data:
        return muts

    for feat in data.get("features", []):
        if feat.get("type") != "Mutagenesis":
            continue
        desc = feat.get("description", "")
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None:
            continue
        if end is None:
            end = start

        alt_seq = feat.get("alternativeSequence", {})
        original = alt_seq.get("originalSequence", "")
        alternatives = alt_seq.get("alternativeSequences", [])
        alt_str = "/".join(alternatives) if alternatives else "?"

        label = f"{original}→{alt_str}: {desc}" if desc else f"{original}→{alt_str}"

        for p in range(int(start), int(end) + 1):
            muts[p].append(label)

    log.info("UniProt mutagenesis: %d positions for %s", len(muts), uniprot)
    return muts


# ── 4. Cancer Gene Census context ───────────────────────────────


def _get_gene_cancer_role(gene_symbol: str, uniprot: str) -> str:
    """
    Determine if the gene is a known oncogene, tumor suppressor, or both.
    Uses UniProt keywords and Open Targets as proxy for Cancer Gene Census
    (which requires COSMIC license for direct access).
    """
    roles = set()

    if uniprot:
        url = f"{UNIPROT_REST}/uniprotkb/{uniprot}.json"
        data = http_get_json(url)
        if data:
            keywords = [kw.get("name", "").lower() for kw in data.get("keywords", [])]
            comments = data.get("comments", [])

            if any("tumor suppressor" in kw or "tumour suppressor" in kw for kw in keywords):
                roles.add("tumor_suppressor")
            if any("oncogene" in kw or "proto-oncogene" in kw for kw in keywords):
                roles.add("oncogene")

            # Check function comments for cancer-related terms
            for comment in comments:
                if comment.get("commentType") == "FUNCTION":
                    texts = comment.get("texts", [])
                    for t in texts:
                        val = t.get("value", "").lower()
                        if "tumor suppressor" in val or "tumour suppressor" in val:
                            roles.add("tumor_suppressor")
                        if "oncogen" in val:
                            roles.add("oncogene")

    return "; ".join(sorted(roles)) if roles else ""


# ── 5. Aggregate cancer tier ────────────────────────────────────


def _classify_cancer_tier(cosmic_count: int, is_hotspot: bool) -> str:
    """
    hotspot  – recurrently mutated (≥5 COSMIC obs or known hotspot)
    moderate – observed somatic mutations (2-4 COSMIC obs)
    rare     – ≤1 observation
    none     – no somatic data
    """
    if is_hotspot or cosmic_count >= 5:
        return "hotspot"
    elif cosmic_count >= 2:
        return "moderate"
    elif cosmic_count >= 1:
        return "rare"
    else:
        return "none"


# =====================================================================
#  PUBLIC PIPELINE
# =====================================================================

COLUMNS = [
    "position",
    "aa",
    "cosmic_ids",
    "cosmic_mutation_count",
    "cosmic_aa_changes",
    "cosmic_cancer_types",
    "is_cancer_hotspot",
    "hotspot_score",
    "cancer_tier",
    "mutagenesis_data",
    "gene_cancer_role",
]


def annotate_cancer(ids: dict, sequence: str) -> list[dict]:
    """Full somatic cancer annotation per residue."""
    seq_len = len(sequence)
    uniprot = ids.get("uniprot")
    ens_prot = ids.get("ensembl_protein")
    gene = ids.get("gene_symbol", "")

    # ── Fetch data ───────────────────────────────────────────
    cosmic = _fetch_cosmic_via_ensembl(ens_prot)
    hotspots = _fetch_cancer_hotspots(gene)
    mutagenesis = _fetch_uniprot_mutagenesis(uniprot) if uniprot else {}
    cancer_role = _get_gene_cancer_role(gene, uniprot)

    # ── Assemble per-residue ─────────────────────────────────
    results = []
    for pos in range(1, seq_len + 1):
        row = {"position": pos, "aa": sequence[pos - 1]}

        cos = cosmic.get(pos, {})
        hs = hotspots.get(pos, {})

        cosmic_ids = cos.get("cosmic_ids", [])
        cosmic_count = cos.get("count", 0)
        aa_changes = cos.get("aa_changes", [])
        cancer_types = cos.get("cancer_types", set()) | hs.get("cancer_types", set())

        is_hotspot = pos in hotspots
        hs_score = hs.get("hotspot_score")

        row["cosmic_ids"] = "; ".join(cosmic_ids) if cosmic_ids else ""
        row["cosmic_mutation_count"] = cosmic_count
        row["cosmic_aa_changes"] = "; ".join(sorted(set(aa_changes))) if aa_changes else ""
        row["cosmic_cancer_types"] = "; ".join(sorted(cancer_types)) if cancer_types else ""
        row["is_cancer_hotspot"] = is_hotspot
        row["hotspot_score"] = hs_score if hs_score is not None else ""
        row["cancer_tier"] = _classify_cancer_tier(cosmic_count, is_hotspot)
        row["mutagenesis_data"] = "; ".join(mutagenesis.get(pos, [])) if pos in mutagenesis else ""
        row["gene_cancer_role"] = cancer_role

        results.append(row)

    n_hotspot = sum(1 for r in results if r["cancer_tier"] == "hotspot")
    n_somatic = sum(1 for r in results if r["cosmic_mutation_count"] > 0)
    log.info(
        "Cancer: %d positions with somatic data, %d hotspots. Gene role: %s",
        n_somatic,
        n_hotspot,
        cancer_role or "none",
    )

    return results


# ── Entry point ──────────────────────────────────────────────────


def main():
    ap = argparse.ArgumentParser(description="Somatic cancer mutation annotations")
    ap.add_argument("--ids-json", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    with open(args.ids_json) as f:
        ids = json.load(f)
    record = SeqIO.read(args.fasta, "fasta")

    rows = annotate_cancer(ids, str(record.seq))

    outpath = Path(args.outdir) / "cancer.tsv"
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    log.info("Wrote %s", outpath)


if __name__ == "__main__":
    main()
