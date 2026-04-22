#!/usr/bin/env python3
"""
open_targets_client.py — Lightweight client for the Open Targets Platform GraphQL API.

Usage:
    python3 open_targets_client.py search-target BRAF
    python3 open_targets_client.py search-disease "breast cancer"
    python3 open_targets_client.py target-info ENSG00000157764
    python3 open_targets_client.py known-drugs ENSG00000157764
    python3 open_targets_client.py association ENSG00000157764 EFO_0000305
    python3 open_targets_client.py validate "BRAF,KRAS" "melanoma"
    python3 open_targets_client.py find-targets "Alzheimer disease" [top_n]
"""

import csv
import json
import os
import sys
from datetime import datetime
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

API_URL = "https://api.platform.opentargets.org/api/v4/graphql"


def graphql_query(query: str, variables: dict | None = None) -> dict:
    payload = json.dumps({"query": query, "variables": variables or {}}).encode()
    req = Request(API_URL, data=payload, headers={"Content-Type": "application/json"})
    try:
        with urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode())
            if "errors" in data:
                print(f"  GraphQL errors: {data['errors']}", file=sys.stderr)
            return data
    except (URLError, HTTPError) as e:
        return {"error": str(e)}


# ── Search ─────────────────────────────────────────────────────────────────


def search_target(name: str) -> dict:
    q = """query($q: String!) {
      search(queryString: $q, entityNames: ["target"], page: {size: 5, index: 0}) {
        hits { id name entity description }
      }
    }"""
    return graphql_query(q, {"q": name})


def search_disease(name: str) -> dict:
    q = """query($q: String!) {
      search(queryString: $q, entityNames: ["disease"], page: {size: 5, index: 0}) {
        hits { id name entity description }
      }
    }"""
    return graphql_query(q, {"q": name})


# ── Target Info (tractability, safety, constraint, pathways) ───────────────


def get_target_info(ensembl_id: str) -> dict:
    q = """query($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        id
        approvedSymbol
        approvedName
        biotype
        functionDescriptions
        subcellularLocations { location }
        tractability { label modality value }
        geneticConstraint { constraintType score oe oeLower oeUpper }
        pathways { pathway }
      }
    }"""
    return graphql_query(q, {"ensemblId": ensembl_id})


# ── Drugs & Clinical Candidates ────────────────────────────────────────────


def get_known_drugs(ensembl_id: str) -> dict:
    q = """query($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        approvedSymbol
        drugAndClinicalCandidates {
          count
          rows {
            id
            maxClinicalStage
            drug { id name drugType }
            diseases { id name }
          }
        }
      }
    }"""
    return graphql_query(q, {"ensemblId": ensembl_id})


# ── Find Targets for Disease ──────────────────────────────────────────────


def find_targets_for_disease(efo_id: str, top_n: int = 10) -> dict:
    q = """query($efoId: String!, $size: Int!) {
      disease(efoId: $efoId) {
        name
        associatedTargets(page: {size: $size, index: 0}) {
          count
          rows {
            target { id approvedSymbol approvedName }
            score
            datatypeScores { id score }
          }
        }
      }
    }"""
    return graphql_query(q, {"efoId": efo_id, "size": top_n})


# ── Association (target-disease) ───────────────────────────────────────────


def get_association(ensembl_id: str, efo_id: str) -> dict:
    q = """query($efoId: String!, $ensemblIds: [String!]) {
      disease(efoId: $efoId) {
        name
        associatedTargets(Bs: $ensemblIds, page: {size: 1, index: 0}) {
          rows {
            target { id approvedSymbol }
            score
            datatypeScores { id score }
            datasourceScores { id score }
          }
        }
      }
    }"""
    result = graphql_query(q, {"efoId": efo_id, "ensemblIds": [ensembl_id]})
    # Check if we got valid rows
    got_rows = False
    if "data" in result and result["data"].get("disease"):
        rows = result["data"]["disease"]["associatedTargets"]["rows"]
        got_rows = len(rows) > 0
    # Fallback: if Bs filter returns nothing or query errored, try without filter
    if not got_rows:
        q2 = """query($efoId: String!) {
          disease(efoId: $efoId) {
            name
            associatedTargets(page: {size: 500, index: 0}) {
              rows {
                target { id approvedSymbol }
                score
                datatypeScores { id score }
                datasourceScores { id score }
              }
            }
          }
        }"""
        result2 = graphql_query(q2, {"efoId": efo_id})
        if "data" in result2 and result2["data"].get("disease"):
            rows = result2["data"]["disease"]["associatedTargets"]["rows"]
            matched = [r for r in rows if r["target"]["id"] == ensembl_id]
            result2["data"]["disease"]["associatedTargets"]["rows"] = matched
            result = result2
    return result


# ── Scoring Logic ──────────────────────────────────────────────────────────


def score_clinical(assoc_data: dict, drugs_data: dict, efo_id: str) -> tuple[int, str]:
    reasons = []

    # Map maxClinicalStage strings to numeric phases
    stage_map = {
        "Phase IV": 4,
        "Approved": 4,
        "Phase III": 3,
        "Phase II": 2,
        "Phase I": 1,
        "Phase I (Early)": 0.5,
    }

    max_phase = 0
    drug_names = []
    dacc = drugs_data.get("data", {}).get("target", {}).get("drugAndClinicalCandidates")
    if dacc:
        for row in dacc.get("rows", []):
            stage = row.get("maxClinicalStage", "")
            ph = stage_map.get(stage, 0)
            drug_obj = row.get("drug") or {}
            drug_name = drug_obj.get("name", "unknown")
            if ph > max_phase:
                max_phase = ph
            if drug_name not in drug_names and len(drug_names) < 5:
                drug_names.append(drug_name)
        n_total = dacc.get("count", len(drug_names))
        if drug_names:
            reasons.append(f"{n_total} drug(s), max phase {max_phase}")

    # Check genetic association from association data
    genetic_score = 0
    if assoc_data.get("data", {}).get("disease"):
        rows = assoc_data["data"]["disease"]["associatedTargets"]["rows"]
        if rows:
            for ds in rows[0].get("datatypeScores", []):
                cid = ds["id"].lower()
                if "genetic" in cid or "somatic" in cid:
                    genetic_score = max(genetic_score, ds["score"])
                    reasons.append(f"Genetic score ({ds['id']}): {ds['score']:.2f}")

    if max_phase >= 4:
        return 5, "; ".join(reasons) or "Approved drug"
    elif max_phase == 3 or genetic_score > 0.5:
        return 4, "; ".join(reasons)
    elif max_phase >= 1 or genetic_score > 0.3:
        return 3, "; ".join(reasons)
    elif genetic_score > 0.1:
        return 2, "; ".join(reasons)
    elif genetic_score > 0 or reasons:
        return 1, "; ".join(reasons)
    return 0, "No clinical or genetic evidence found"


def score_druggability(target_data: dict, drugs_data: dict) -> tuple[int, str]:
    reasons = []
    tract_sm = False
    tract_ab = False

    sm_labels = []
    ab_labels = []
    if target_data.get("data", {}).get("target"):
        t = target_data["data"]["target"]
        for tr in t.get("tractability", []):
            if tr.get("value"):
                if tr["modality"] == "SM":
                    tract_sm = True
                    sm_labels.append(tr["label"])
                elif tr["modality"] == "AB":
                    tract_ab = True
                    ab_labels.append(tr["label"])
    if sm_labels:
        reasons.append(f"SM({len(sm_labels)}): {', '.join(sm_labels[:2])}")
    if ab_labels:
        reasons.append(f"AB({len(ab_labels)}): {', '.join(ab_labels[:2])}")

    stage_map = {
        "Phase IV": 4,
        "Approved": 4,
        "Phase III": 3,
        "Phase II": 2,
        "Phase I": 1,
        "Phase I (Early)": 0.5,
    }
    has_approved = False
    has_clinical = False
    n_drugs = 0
    dacc = drugs_data.get("data", {}).get("target", {}).get("drugAndClinicalCandidates")
    if dacc:
        n_drugs = dacc.get("count", 0)
        for row in dacc.get("rows", []):
            stage = row.get("maxClinicalStage", "")
            ph = stage_map.get(stage, 0)
            if ph >= 4:
                has_approved = True
            elif ph >= 1:
                has_clinical = True
        if n_drugs > 0:
            reasons.append(f"{n_drugs} drug/clinical candidate(s)")

    if has_approved:
        return 5, "; ".join(reasons)
    elif has_clinical and (tract_sm or tract_ab):
        return 4, "; ".join(reasons)
    elif tract_sm or tract_ab:
        return 3, "; ".join(reasons)
    elif has_clinical:
        return 3, "; ".join(reasons)
    elif n_drugs > 0:
        return 2, "; ".join(reasons)
    return 0, "No tractability or drug data"


def score_pathway(target_data: dict, assoc_data: dict) -> tuple[int, str]:
    reasons = []
    pathway_count = 0

    if target_data.get("data", {}).get("target"):
        t = target_data["data"]["target"]
        pathways = t.get("pathways") or []
        pathway_count = len(pathways)
        if pathway_count > 0:
            top = [p.get("pathway", "")[:40] for p in pathways[:2]]
            reasons.append(f"{pathway_count} pw: {'; '.join(top)}")

    lit_score = 0
    expr_score = 0
    if assoc_data.get("data", {}).get("disease"):
        rows = assoc_data["data"]["disease"]["associatedTargets"]["rows"]
        if rows:
            for ds in rows[0].get("datatypeScores", []):
                cid = ds["id"].lower()
                if "literature" in cid:
                    lit_score = ds["score"]
                    reasons.append(f"Literature: {ds['score']:.2f}")
                if "expression" in cid or "rna" in cid:
                    expr_score = ds["score"]
                    reasons.append(f"Expression: {ds['score']:.2f}")

    if pathway_count >= 5 and (expr_score > 0.1 or lit_score > 0.3):
        return 5, "; ".join(reasons)
    elif pathway_count >= 3 or (pathway_count >= 1 and lit_score > 0.2):
        return 4, "; ".join(reasons)
    elif pathway_count >= 1 and (expr_score > 0 or lit_score > 0):
        return 3, "; ".join(reasons)
    elif pathway_count >= 1 or lit_score > 0 or expr_score > 0:
        return 2, "; ".join(reasons)
    elif reasons:
        return 1, "; ".join(reasons)
    return 0, "No pathway or expression data"


def composite_score(clinical: int, druggability: int, pathway: int) -> float:
    return round(clinical * 0.40 + druggability * 0.35 + pathway * 0.25, 2)


def confidence_level(comp: float, n_sources: int) -> str:
    if comp >= 3.5 and n_sources >= 3:
        return "High"
    elif comp >= 2.0:
        return "Medium"
    return "Low"


# ── Validation Pipeline ───────────────────────────────────────────────────


def validate(targets: list[str], diseases: list[str]) -> dict:
    results = []

    # Resolve targets
    resolved_targets = {}
    for t in targets:
        if t.startswith("ENSG"):
            resolved_targets[t] = t
        else:
            r = search_target(t)
            hits = r.get("data", {}).get("search", {}).get("hits", [])
            if hits:
                resolved_targets[t] = hits[0]["id"]
                print(f"  {t} → {hits[0]['id']}")
            else:
                print(f"  WARN: target '{t}' not found")

    # Resolve diseases
    resolved_diseases = {}
    for d in diseases:
        if any(d.startswith(p) for p in ["EFO_", "MONDO_", "HP_", "Orphanet_"]):
            resolved_diseases[d] = d
        else:
            r = search_disease(d)
            hits = r.get("data", {}).get("search", {}).get("hits", [])
            if hits:
                resolved_diseases[d] = hits[0]["id"]
                print(f"  {d} → {hits[0]['id']}")
            else:
                print(f"  WARN: disease '{d}' not found")

    raw_data = {}
    for t_name, t_id in resolved_targets.items():
        for d_name, d_id in resolved_diseases.items():
            pair_key = f"{t_name}|{d_name}"
            print(f"  scoring {t_name}×{d_name}")

            assoc = get_association(t_id, d_id)
            target_info = get_target_info(t_id)
            drugs = get_known_drugs(t_id)

            raw_data[pair_key] = {
                "association": assoc,
                "target_info": target_info,
                "known_drugs": drugs,
            }

            clin_score, clin_reason = score_clinical(assoc, drugs, d_id)
            drug_score, drug_reason = score_druggability(target_info, drugs)
            path_score, path_reason = score_pathway(target_info, assoc)
            comp = composite_score(clin_score, drug_score, path_score)
            n_sources = sum(1 for s in [clin_score, drug_score, path_score] if s > 0)
            conf = confidence_level(comp, n_sources)

            symbol = t_name
            full_name = ""
            context = ""
            if target_info.get("data", {}).get("target"):
                ti = target_info["data"]["target"]
                symbol = ti.get("approvedSymbol", t_name)
                full_name = ti.get("approvedName", "")
                funcs = ti.get("functionDescriptions") or []
                context = funcs[0][:80] if funcs else ""

            disease_name = d_name
            if assoc.get("data", {}).get("disease"):
                disease_name = assoc["data"]["disease"].get("name", d_name)

            # Safety score
            safety = 5
            safety_reason = "No safety liabilities reported"
            if target_info.get("data", {}).get("target"):
                ti = target_info["data"]["target"]
                gc = ti.get("geneticConstraint") or []
                for g in gc:
                    if g.get("constraintType") == "lof" and g.get("score") is not None:
                        if g["score"] > 0.9:
                            safety = 1
                            safety_reason = (
                                f"High LoF intolerance (pLI-like score={g['score']:.2f})"
                            )

            results.append(
                {
                    "target": symbol,
                    "target_id": t_id,
                    "full_name": full_name,
                    "context": context,
                    "disease": disease_name,
                    "disease_id": d_id,
                    "clinical": clin_score,
                    "clinical_reason": clin_reason,
                    "druggability": drug_score,
                    "druggability_reason": drug_reason,
                    "pathway": path_score,
                    "pathway_reason": path_reason,
                    "safety": safety,
                    "safety_reason": safety_reason,
                    "composite": comp,
                    "confidence": conf,
                }
            )

    results.sort(key=lambda x: x["composite"], reverse=True)
    return {"results": results, "raw_data": raw_data}


# ── Output Formatters ──────────────────────────────────────────────────────


def format_markdown_table(results):
    lines = [
        "| Rank | Target | Disease | Clinical | Druggability | Pathway | Composite | Confidence |",
        "|------|--------|---------|:--------:|:----------:|:-------:|:---------:|:----------:|",
    ]
    for i, r in enumerate(results, 1):
        lines.append(
            f"| {i} | **{r['target']}** | {r['disease']} | {r['clinical']} | "
            f"{r['druggability']} | {r['pathway']} | **{r['composite']}** | {r['confidence']} |"
        )
    return "\n".join(lines)


def format_narrative(results):
    sections = []
    for r in results:
        sections.append(
            f"**{r['target']}×{r['disease']}** C{r['clinical']} D{r['druggability']} P{r['pathway']} S{r['safety']} → {r['composite']} ({r['confidence']})\n"
            f"  Clin: {r['clinical_reason'][:120]}\n"
            f"  Drug: {r['druggability_reason'][:120]}\n"
            f"  Path: {r['pathway_reason'][:120]}\n"
            f"  Safe: {r['safety_reason'][:80]}"
        )
    return "\n".join(sections)


def generate_html(results, disease_name):
    targets_for_js = []
    for r in results:
        targets_for_js.append(
            {
                "name": r["target"],
                "fullName": r.get("full_name", ""),
                "ensemblId": r.get("target_id", ""),
                "disease": r.get("disease", ""),
                "context": r.get("context", ""),
                "scores": {
                    "genetic": min(r["clinical"], 5),
                    "trials": min(r["clinical"], 5),
                    "clinical_total": r["clinical"],
                    "sm_tractability": min(r["druggability"], 5),
                    "ab_tractability": max(0, r["druggability"] - 1),
                    "druggability_total": r["druggability"],
                    "pathways": min(r["pathway"], 5),
                    "expression": max(0, r["pathway"] - 1),
                    "literature": min(r["pathway"], 5),
                    "pathway_total": r["pathway"],
                    "safety": r.get("safety", 5),
                },
                "reasons": {
                    "genetic": r.get("clinical_reason", ""),
                    "trials": r.get("clinical_reason", ""),
                    "clinical_total": r.get("clinical_reason", ""),
                    "sm_tractability": r.get("druggability_reason", ""),
                    "ab_tractability": r.get("druggability_reason", ""),
                    "druggability_total": r.get("druggability_reason", ""),
                    "pathways": r.get("pathway_reason", ""),
                    "expression": r.get("pathway_reason", ""),
                    "literature": r.get("pathway_reason", ""),
                    "pathway_total": r.get("pathway_reason", ""),
                    "safety": r.get("safety_reason", "No safety liabilities reported"),
                },
                "notes": (
                    f"Clinical: {r.get('clinical_reason', 'N/A')}. "
                    f"Druggability: {r.get('druggability_reason', 'N/A')}. "
                    f"Pathway: {r.get('pathway_reason', 'N/A')}."
                ),
            }
        )

    targets_json = json.dumps(targets_for_js, indent=2)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    unique_targets = list(dict.fromkeys(r["target"] for r in results))
    unique_diseases = list(dict.fromkeys(r["disease"] for r in results))
    subtitle = (
        f"{len(unique_targets)} target(s) ({', '.join(unique_targets)}) scored against "
        f"{len(unique_diseases)} disease(s) ({', '.join(unique_diseases)}) across 7 evidence sub-dimensions. "
        f"Scores 0-5 per dimension. Weighted composite determines ranking. "
        f"Safety score of 0 = veto. Hover scores for details. Click target names for notes."
    )

    template_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "report_template.html")
    if not os.path.exists(template_path):
        template_path = os.path.join(".", "report_template.html")
    if not os.path.exists(template_path):
        return "<!-- ERROR: report_template.html not found -->"

    with open(template_path, "r") as f:
        html = f.read()

    html = html.replace("{{TARGETS_JSON}}", targets_json)
    html = html.replace("{{DISEASE_NAME}}", disease_name)
    html = html.replace("{{SUBTITLE}}", subtitle)
    html = html.replace("{{TIMESTAMP}}", timestamp)
    return html


def generate_pathway_html(all_pathway_data: dict) -> str:
    """Generate an interactive HTML pathway visualization for multiple targets."""
    data_json = json.dumps(all_pathway_data, indent=2)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    targets_list = list(all_pathway_data.keys())

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Pathway Map — {", ".join(targets_list)}</title>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;600;700&display=swap" rel="stylesheet">
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:'IBM Plex Sans',system-ui,sans-serif;background:#f8f9fb;color:#374151;padding:24px}}
.container{{max-width:1100px;margin:0 auto;background:white;border-radius:12px;padding:28px;box-shadow:0 2px 12px rgba(0,0,0,.06)}}
h1{{font-size:20px;font-weight:700;color:#1e3a5f;margin-bottom:4px}}
.sub{{font-size:13px;color:#6b7280;margin-bottom:20px}}
.controls{{display:flex;gap:10px;flex-wrap:wrap;margin-bottom:20px}}
.chip{{padding:6px 16px;border-radius:20px;font-size:13px;font-weight:600;cursor:pointer;border:1.5px solid #e5e7eb;background:white;transition:all .2s}}
.chip.active{{background:#1e3a5f;color:white;border-color:#1e3a5f}}
.chip:hover{{border-color:#1e3a5f}}
.filter-bar{{display:flex;gap:10px;align-items:center;margin-bottom:16px}}
.filter-bar input{{flex:1;padding:8px 14px;border:1px solid #e5e7eb;border-radius:8px;font-size:13px;outline:none}}
.filter-bar input:focus{{border-color:#3b82f6;box-shadow:0 0 0 3px rgba(59,130,246,.1)}}
.stats{{font-size:12px;color:#9ca3af;margin-bottom:16px}}
svg{{width:100%;display:block}}
.pw-node{{cursor:pointer;transition:opacity .15s}}
.pw-node:hover{{opacity:.85}}
.pw-label{{font-family:'IBM Plex Sans',sans-serif;fill:#374151;font-size:12px}}
.pw-label-bold{{font-family:'IBM Plex Mono',monospace;fill:#1e3a5f;font-size:13px;font-weight:700}}
.target-circle{{stroke-width:2.5}}
.pathway-rect{{stroke-width:.5;rx:6}}
.edge{{stroke-width:1;fill:none;opacity:.4}}
.edge.highlight{{opacity:.9;stroke-width:2}}
.legend{{margin-top:20px;padding:16px;background:#f9fafb;border-radius:8px;border:1px solid #e5e7eb;font-size:12px;color:#6b7280}}
.legend-items{{display:flex;gap:16px;flex-wrap:wrap;margin-top:8px}}
.legend-item{{display:flex;align-items:center;gap:6px}}
.legend-dot{{width:14px;height:14px;border-radius:50%;border:2px solid}}
.detail-panel{{margin-top:16px;padding:16px;background:#f0f5ff;border-radius:8px;border-left:4px solid #1e3a5f;display:none}}
.detail-panel.visible{{display:block}}
.detail-title{{font-family:'IBM Plex Mono',monospace;font-weight:700;color:#1e3a5f;margin-bottom:6px}}
.detail-body{{font-size:13px;color:#374151;line-height:1.6}}
.detail-body a{{color:#3b82f6;text-decoration:none}}
.detail-body a:hover{{text-decoration:underline}}
.shared-badge{{display:inline-block;padding:2px 8px;border-radius:10px;font-size:10px;font-weight:700;background:#fef3c7;color:#92400e;margin-left:6px}}
.generated{{margin-top:16px;font-size:11px;color:#9ca3af;text-align:right}}
</style>
</head>
<body>
<div class="container">
  <h1>Pathway Map</h1>
  <p class="sub">{len(targets_list)} target(s): {", ".join(targets_list)}. Pathways from Reactome via Open Targets Platform. Click a pathway for details. Shared pathways are highlighted.</p>

  <div class="controls" id="controls"></div>
  <div class="filter-bar">
    <input type="text" id="search" placeholder="Filter pathways..." oninput="filterPathways(this.value)">
  </div>
  <div class="stats" id="stats"></div>

  <div id="viz"></div>

  <div class="detail-panel" id="detail-panel">
    <div class="detail-title" id="detail-title"></div>
    <div class="detail-body" id="detail-body"></div>
  </div>

  <div class="legend">
    <strong>Legend</strong>
    <div class="legend-items" id="legend-items"></div>
    <div style="margin-top:8px">Pathways appearing in multiple targets are marked <span class="shared-badge">SHARED</span>. Click any pathway node to see details and a Reactome link.</div>
  </div>
  <div class="generated">Generated {timestamp} · Data: Open Targets Platform + Reactome</div>
</div>

<script>
const DATA = {data_json};

const COLORS = [
  {{fill:'#dbeafe',stroke:'#3b82f6',text:'#1e40af',name:'Blue'}},
  {{fill:'#d1fae5',stroke:'#10b981',text:'#065f46',name:'Green'}},
  {{fill:'#fef3c7',stroke:'#f59e0b',text:'#92400e',name:'Amber'}},
  {{fill:'#fce7f3',stroke:'#ec4899',text:'#9d174d',name:'Pink'}},
  {{fill:'#e0e7ff',stroke:'#6366f1',text:'#3730a3',name:'Indigo'}},
  {{fill:'#fde8e8',stroke:'#ef4444',text:'#991b1b',name:'Red'}},
];

const targets = Object.keys(DATA);
let activeTargets = new Set(targets);
let filterText = '';
let selectedPathway = null;

function getColor(i) {{ return COLORS[i % COLORS.length]; }}

function buildPathwayIndex() {{
  const idx = {{}};
  targets.forEach((t, ti) => {{
    (DATA[t].pathways || []).forEach(p => {{
      const key = p.id || p.term;
      if (!idx[key]) idx[key] = {{id: p.id, term: p.term, targets: []}};
      if (!idx[key].targets.includes(t)) idx[key].targets.push(t);
    }});
  }});
  return idx;
}}

function render() {{
  const idx = buildPathwayIndex();
  let pathways = Object.values(idx).filter(p =>
    p.targets.some(t => activeTargets.has(t)) &&
    (!filterText || p.term.toLowerCase().includes(filterText.toLowerCase()))
  );
  pathways.sort((a, b) => b.targets.length - a.targets.length || a.term.localeCompare(b.term));

  const shared = pathways.filter(p => p.targets.length > 1).length;
  document.getElementById('stats').textContent =
    `${{pathways.length}} pathways shown (${{shared}} shared across targets)`;

  // Controls
  const cDiv = document.getElementById('controls');
  cDiv.innerHTML = '';
  targets.forEach((t, i) => {{
    const c = document.createElement('span');
    c.className = 'chip' + (activeTargets.has(t) ? ' active' : '');
    c.style.borderColor = activeTargets.has(t) ? getColor(i).stroke : '';
    c.style.background = activeTargets.has(t) ? getColor(i).fill : '';
    c.style.color = activeTargets.has(t) ? getColor(i).text : '';
    c.textContent = t;
    c.onclick = () => {{
      if (activeTargets.has(t)) activeTargets.delete(t); else activeTargets.add(t);
      render();
    }};
    cDiv.appendChild(c);
  }});

  // Legend
  const leg = document.getElementById('legend-items');
  leg.innerHTML = '';
  targets.forEach((t, i) => {{
    const d = document.createElement('div');
    d.className = 'legend-item';
    d.innerHTML = `<div class="legend-dot" style="background:${{getColor(i).fill}};border-color:${{getColor(i).stroke}}"></div>${{t}}`;
    leg.appendChild(d);
  }});

  // SVG
  const targetSpacing = 180;
  const totalTargetWidth = targets.length * targetSpacing;
  const svgW = Math.max(800, totalTargetWidth + 100);
  const rowH = 28;
  const gapH = 4;
  const topMargin = 80;
  const pathwayStartY = topMargin + 60;
  const svgH = pathwayStartY + pathways.length * (rowH + gapH) + 40;
  const targetY = 40;
  const pwX = 220;
  const pwW = svgW - pwX - 20;

  let svg = `<svg viewBox="0 0 ${{svgW}} ${{svgH}}" xmlns="http://www.w3.org/2000/svg">`;

  // Target nodes (circles at top)
  const targetPositions = {{}};
  targets.forEach((t, i) => {{
    const x = pwX + (i + 0.5) * (pwW / targets.length);
    targetPositions[t] = {{x, y: targetY}};
    const col = getColor(i);
    const opacity = activeTargets.has(t) ? 1 : 0.3;
    svg += `<circle cx="${{x}}" cy="${{targetY}}" r="22" fill="${{col.fill}}" stroke="${{col.stroke}}" class="target-circle" opacity="${{opacity}}"/>`;
    svg += `<text x="${{x}}" y="${{targetY + 1}}" text-anchor="middle" dominant-baseline="central" class="pw-label-bold" fill="${{col.text}}" opacity="${{opacity}}">${{t}}</text>`;
  }});

  // Pathway rows
  pathways.forEach((p, pi) => {{
    const y = pathwayStartY + pi * (rowH + gapH);
    const isShared = p.targets.length > 1;
    const isSelected = selectedPathway === p.id;
    const bgFill = isSelected ? '#e0e7ff' : (isShared ? '#fefce8' : '#fafbfc');
    const borderCol = isSelected ? '#6366f1' : (isShared ? '#fde047' : '#e5e7eb');

    svg += `<g class="pw-node" onclick="selectPathway('${{p.id.replace(/'/g, "\\\\'")}}')">`;
    svg += `<rect x="4" y="${{y}}" width="${{svgW - 8}}" height="${{rowH}}" fill="${{bgFill}}" stroke="${{borderCol}}" class="pathway-rect"/>`;
    const label = p.term.length > 45 ? p.term.slice(0, 42) + '...' : p.term;
    svg += `<text x="14" y="${{y + rowH / 2}}" dominant-baseline="central" class="pw-label">${{label}}</text>`;
    if (isShared) svg += `<text x="${{pwX - 10}}" y="${{y + rowH / 2}}" dominant-baseline="central" text-anchor="end" font-size="9" font-weight="700" fill="#92400e">${{p.targets.length}}</text>`;

    // Edges from target to this pathway
    p.targets.forEach(t => {{
      if (!activeTargets.has(t)) return;
      const ti = targets.indexOf(t);
      const col = getColor(ti);
      const tx = targetPositions[t].x;
      const py = y + rowH / 2;
      const midY = topMargin + 30;
      const highlight = isSelected ? ' highlight' : '';
      svg += `<path d="M${{tx}} ${{targetY + 24}} Q${{tx}} ${{midY}} ${{tx}} ${{midY}} L${{tx}} ${{py}}" stroke="${{col.stroke}}" class="edge${{highlight}}" marker-end="none"/>`;
      // Dot on the pathway row
      svg += `<circle cx="${{tx}}" cy="${{py}}" r="4" fill="${{col.stroke}}" opacity="0.7"/>`;
    }});

    svg += `</g>`;
  }});

  svg += '</svg>';
  document.getElementById('viz').innerHTML = svg;

  // Detail panel
  const panel = document.getElementById('detail-panel');
  if (selectedPathway) {{
    const p = Object.values(idx).find(x => x.id === selectedPathway);
    if (p) {{
      document.getElementById('detail-title').textContent = p.term;
      let html = `<strong>Reactome ID:</strong> ${{p.id}}<br>`;
      html += `<strong>Targets:</strong> ${{p.targets.join(', ')}}<br>`;
      html += `<a href="https://reactome.org/content/detail/${{p.id}}" target="_blank">View in Reactome &rarr;</a>`;
      document.getElementById('detail-body').innerHTML = html;
      panel.classList.add('visible');
    }}
  }} else {{
    panel.classList.remove('visible');
  }}
}}

function selectPathway(id) {{
  selectedPathway = selectedPathway === id ? null : id;
  render();
}}

function filterPathways(text) {{
  filterText = text;
  render();
}}

render();
</script>
</body>
</html>"""


def _make_slug(results):
    targets = list(dict.fromkeys(r["target"] for r in results))
    diseases = list(dict.fromkeys(r["disease"] for r in results))

    def slugify(s):
        return s.lower().replace(" ", "_").replace("'", "")

    return f"{'_'.join(slugify(t) for t in targets)}_vs_{'_'.join(slugify(d) for d in diseases)}"


def save_results(validation, output_dir="results"):
    os.makedirs(output_dir, exist_ok=True)
    results = validation["results"]
    slug = _make_slug(results) if results else "unknown"

    with open(f"{output_dir}/validation_{slug}.md", "w") as f:
        f.write(
            f"# Target Validation Report\n**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n"
        )
        f.write(
            f"## Summary\n\n{format_markdown_table(results)}\n\n## Details\n\n{format_narrative(results)}\n"
        )

    if results:
        fieldnames = [k for k in results[0].keys() if k != "full_name" and k != "context"]
        with open(f"{output_dir}/scores_{slug}.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
            w.writeheader()
            w.writerows(results)

    with open(f"{output_dir}/raw_data_{slug}.json", "w") as f:
        json.dump(validation["raw_data"], f, indent=2, default=str)

    disease_label = (
        ", ".join(dict.fromkeys(r["disease"] for r in results)) if results else "unknown"
    )
    html = generate_html(results, disease_label)
    html_path = f"{output_dir}/validation_{slug}.html"
    with open(html_path, "w") as f:
        f.write(html)

    print(f"\n  Files saved to {output_dir}/")
    print(f"  Open {html_path} in a browser to view the interactive matrix.")


# ── CLI ────────────────────────────────────────────────────────────────────


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)

    cmd = sys.argv[1]

    if cmd == "search-target" and len(sys.argv) >= 3:
        print(json.dumps(search_target(sys.argv[2]), indent=2))
    elif cmd == "search-disease" and len(sys.argv) >= 3:
        print(json.dumps(search_disease(" ".join(sys.argv[2:])), indent=2))
    elif cmd == "association" and len(sys.argv) >= 4:
        print(json.dumps(get_association(sys.argv[2], sys.argv[3]), indent=2))
    elif cmd == "target-info" and len(sys.argv) >= 3:
        print(json.dumps(get_target_info(sys.argv[2]), indent=2))
    elif cmd == "known-drugs" and len(sys.argv) >= 3:
        print(json.dumps(get_known_drugs(sys.argv[2]), indent=2))
    elif cmd == "find-targets" and len(sys.argv) >= 3:
        disease_query = " ".join(sys.argv[2:])
        top_n = 10
        # Check if last arg is a number (top_n)
        if sys.argv[-1].isdigit():
            top_n = int(sys.argv[-1])
            disease_query = " ".join(sys.argv[2:-1])
        # Resolve disease
        efo_id = disease_query
        if not any(disease_query.startswith(p) for p in ["EFO_", "MONDO_", "HP_", "Orphanet_"]):
            r = search_disease(disease_query)
            hits = r.get("data", {}).get("search", {}).get("hits", [])
            if hits:
                efo_id = hits[0]["id"]
                print(f"  {disease_query} → {efo_id}")
            else:
                print(f"  WARN: disease '{disease_query}' not found")
                sys.exit(1)
        result = find_targets_for_disease(efo_id, top_n)
        disease_data = result.get("data", {}).get("disease")
        if not disease_data:
            print("  No data returned for this disease.")
            sys.exit(1)
        rows = disease_data["associatedTargets"]["rows"]
        total = disease_data["associatedTargets"]["count"]
        print(f"  {disease_data['name']}: {total} associated targets (showing top {len(rows)})\n")
        print("| Rank | Symbol | Name | Score | Top datatype |")
        print("|------|--------|------|:-----:|--------------|")
        for i, row in enumerate(rows, 1):
            t = row["target"]
            score = row["score"]
            dts = sorted(row.get("datatypeScores", []), key=lambda d: d["score"], reverse=True)
            top_dt = f"{dts[0]['id']} ({dts[0]['score']:.2f})" if dts else "—"
            name = (t.get("approvedName") or "")[:40]
            print(f"| {i} | **{t['approvedSymbol']}** | {name} | {score:.3f} | {top_dt} |")
    elif cmd == "validate" and len(sys.argv) >= 4:
        targets = [t.strip() for t in sys.argv[2].split(",")]
        diseases = [d.strip() for d in sys.argv[3].split(",")]
        print(f"Validating {len(targets)} target(s) × {len(diseases)} disease(s)...\n")
        v = validate(targets, diseases)
        print("\n" + format_markdown_table(v["results"]))
        print("\n" + format_narrative(v["results"]))
        save_results(v)
    elif cmd == "pathways" and len(sys.argv) >= 3:
        targets = [t.strip() for t in sys.argv[2].split(",")]
        print(f"Fetching pathways for {len(targets)} target(s)...\n")
        all_pathway_data = {}
        for t in targets:
            ensembl_id = t
            symbol = t
            if not t.startswith("ENSG"):
                r = search_target(t)
                hits = r.get("data", {}).get("search", {}).get("hits", [])
                if hits:
                    ensembl_id = hits[0]["id"]
                    symbol = hits[0].get("name", t)
                    print(f"  Resolved {t} → {ensembl_id} ({symbol})")
                else:
                    print(f"  WARNING: Could not resolve '{t}'")
                    continue
            info = get_target_info(ensembl_id)
            target_obj = info.get("data", {}).get("target", {})
            pathways = target_obj.get("pathways") or []
            approved = target_obj.get("approvedSymbol", symbol)
            full_name = target_obj.get("approvedName", "")
            all_pathway_data[approved] = {
                "ensemblId": ensembl_id,
                "symbol": approved,
                "fullName": full_name,
                "pathways": [
                    {"id": p.get("pathway", ""), "term": p.get("pathway", "")} for p in pathways
                ],
            }
            print(f"  {approved}: {len(pathways)} pathway(s)")
        # Generate HTML
        os.makedirs("results", exist_ok=True)
        html = generate_pathway_html(all_pathway_data)
        out_path = "results/pathways.html"
        with open(out_path, "w") as f:
            f.write(html)
        print(f"\n  Saved to {out_path}")
        for sym, data in all_pathway_data.items():
            print(f"  {sym}: {len(data['pathways'])} pathways")
    else:
        print(__doc__)
        sys.exit(1)


if __name__ == "__main__":
    main()
