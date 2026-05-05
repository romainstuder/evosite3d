#!/usr/bin/env python3
"""
modules/disease.py
===================
Per-residue genetic disease associations:

  1. UniProt variants   – "Natural variant" features with disease links
  2. ClinVar/Ensembl    – pathogenic missense variants via Ensembl overlap
  3. Open Targets       – gene-level disease associations (context)

Inputs:  resolved_ids.json, sequence.fasta
Outputs: disease.tsv  (position, disease_association, gene_diseases)
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


# ── 1. UniProt disease variants ──────────────────────────────────


def fetch_uniprot_disease_variants(uniprot: str) -> dict[int, list[str]]:
    variants: dict[int, list[str]] = defaultdict(list)
    url = f"{UNIPROT_REST}/uniprotkb/{uniprot}.json"
    data = http_get_json(url)
    if not data:
        return variants

    _DISEASE_KW = {
        "cancer",
        "carcinoma",
        "tumor",
        "tumour",
        "leukemia",
        "deficiency",
        "disorder",
        "dystrophy",
        "anemia",
        "susceptibility",
        "pathogenic",
        "syndrome",
        "disease",
    }

    for feat in data.get("features", []):
        if feat.get("type") != "Natural variant":
            continue

        desc = feat.get("description", "")
        loc = feat.get("location", {})
        pos = loc.get("start", {}).get("value")
        if pos is None:
            continue
        pos = int(pos)

        alt_seq = feat.get("alternativeSequence", {})
        original = alt_seq.get("originalSequence", "?")
        alternatives = alt_seq.get("alternativeSequences", [])
        alt_str = "/".join(alternatives) if alternatives else "?"

        # Detect disease association
        desc_lower = desc.lower()
        is_disease = False
        disease_name = ""

        if " in " in desc_lower:
            is_disease = True
            disease_name = desc[desc_lower.index(" in ") + 4 :].strip().rstrip(".")
        elif any(kw in desc_lower for kw in _DISEASE_KW):
            is_disease = True
            disease_name = desc.strip()

        # Check evidence sources
        for ev in feat.get("evidences", []):
            source = ev.get("source", {})
            if isinstance(source, str):
                source = {"name": source}
            src_name = source.get("name", "").lower()
            if src_name in ("omim", "clinvar"):
                is_disease = True
                if not disease_name:
                    disease_name = source.get("id", "")

        if is_disease:
            label = f"{original}{pos}{alt_str}"
            if disease_name:
                label += f": {disease_name}"
            label += " (UniProt)"
            variants[pos].append(label)

    log.info("UniProt disease variants: %d positions for %s", len(variants), uniprot)
    return variants


# ── 2. ClinVar via Ensembl ───────────────────────────────────────


def fetch_clinvar_variants(ensembl_protein: str) -> dict[int, list[str]]:
    variants: dict[int, list[str]] = defaultdict(list)
    if not ensembl_protein:
        return variants

    data = fetch_translation_overlap(ensembl_protein, feature="transcript_variation")
    if not data:
        return variants

    _PATHOGENIC = {
        "pathogenic",
        "likely_pathogenic",
        "pathogenic/likely_pathogenic",
        "likely pathogenic",
        "risk_factor",
        "association",
    }

    for var in data:
        if "missense_variant" not in var.get("consequence_terms", []):
            continue
        clinical = var.get("clinical_significance", [])
        if not any(c.lower() in _PATHOGENIC for c in clinical):
            continue

        start = var.get("start")
        if start is None:
            continue

        residues = var.get("residues", "")
        var_id = var.get("id", "")
        clin_str = "/".join(clinical)
        label = f"{residues} {var_id} [{clin_str}] (ClinVar/Ensembl)"
        variants[int(start)].append(label)

    log.info("ClinVar variants: %d positions", len(variants))
    return variants


# ── 3. Gene-level disease (Open Targets) ─────────────────────────


def fetch_gene_diseases(ensembl_gene: str) -> list[str]:
    if not ensembl_gene:
        return []
    # url = (
    #     f"https://api.platform.opentargets.org/api/v4/graphql"
    #     f'?query={{target(ensemblId:"{ensembl_gene}"){{associatedDiseases(page:{{size:10}}){{rows{{disease{{name}}score}}}}}}}}'
    # )
    # Simpler: try the REST-like endpoint
    url_rest = (
        f"https://api.platform.opentargets.org/api/v4/target/{ensembl_gene}/associations?size=10"
    )
    data = http_get_json(url_rest)
    diseases = []
    if data and "data" in data:
        for assoc in data["data"]:
            name = assoc.get("disease", {}).get("name", "")
            score = assoc.get("score", 0)
            if name and score > 0.3:
                diseases.append(f"{name} (score={score:.2f})")
    return diseases


# ── Pipeline ─────────────────────────────────────────────────────


def annotate_disease(ids: dict, sequence: str) -> list[dict]:
    seq_len = len(sequence)
    results = [
        {"position": i + 1, "aa": sequence[i], "disease_association": "", "gene_diseases": ""}
        for i in range(seq_len)
    ]

    disease_map: dict[int, list[str]] = defaultdict(list)

    if ids.get("uniprot"):
        up = fetch_uniprot_disease_variants(ids["uniprot"])
        for pos, labels in up.items():
            disease_map[pos].extend(labels)

    if ids.get("ensembl_protein"):
        cv = fetch_clinvar_variants(ids["ensembl_protein"])
        for pos, labels in cv.items():
            disease_map[pos].extend(labels)

    gene_diseases_str = ""
    if ids.get("ensembl_gene"):
        gd = fetch_gene_diseases(ids["ensembl_gene"])
        if gd:
            gene_diseases_str = "; ".join(gd[:5])

    for r in results:
        pos = r["position"]
        if pos in disease_map:
            r["disease_association"] = "; ".join(disease_map[pos])
        r["gene_diseases"] = gene_diseases_str

    n = sum(1 for r in results if r["disease_association"])
    log.info("Disease: %d / %d residues with variant-level annotations.", n, seq_len)
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

    rows = annotate_disease(ids, str(record.seq))

    outpath = Path(args.outdir) / "disease.tsv"
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["position", "aa", "disease_association", "gene_diseases"], delimiter="\t"
        )
        w.writeheader()
        w.writerows(rows)
    log.info("Wrote %s", outpath)


if __name__ == "__main__":
    main()
