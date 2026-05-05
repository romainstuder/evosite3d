#!/usr/bin/env python3
"""
bin/batch_runner.py
====================
Batch mode: run the characterization pipeline on multiple proteins
and generate cross-protein comparative analysis.

  ┌──────────────────────────────────────────────────────────────┐
  │  BATCH MODE                                                  │
  │                                                              │
  │  1. Input: gene list (file or CLI), protein family, or      │
  │     pathway from Reactome/KEGG                              │
  │  2. Runs Nextflow pipeline for each protein                 │
  │  3. Cross-protein comparison:                               │
  │     - Designability index distribution per protein          │
  │     - Shared domains and their conservation patterns        │
  │     - Interface residue overlap between family members      │
  │     - Cancer hotspot enrichment comparison                  │
  │     - Best engineering candidates across the family         │
  │  4. Outputs:                                                 │
  │     - Per-protein TSVs (as usual)                           │
  │     - batch_comparison.tsv (cross-protein summary)          │
  │     - family_report.tsv (per-residue alignment of DI)       │
  │     - engineering_candidates.tsv (ranked mutable sites)      │
  └──────────────────────────────────────────────────────────────┘

Usage:
    # From gene list
    python -m protein_characteriser.batch --genes TP53 BRCA1 EGFR KRAS

    # From file (one gene per line)
    python -m protein_characteriser.batch --gene-list genes.txt

    # From pathway
    python -m protein_characteriser.batch --pathway "p53 signaling"

    # Compare existing results
    python -m protein_characteriser.batch --results-dir results/
"""

import argparse
import csv
import logging
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from api_utils import http_get_json

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


# ── Pathway/family gene lookup ───────────────────────────────────


def fetch_pathway_genes(pathway_query: str) -> list[str]:
    """
    Fetch gene list from Reactome by pathway name search.
    Returns list of HGNC gene symbols.
    """
    # Reactome search
    url = f"https://reactome.org/ContentService/search/query?query={pathway_query}&species=Homo+sapiens&types=Pathway&cluster=true"
    data = http_get_json(url)

    if not data or "results" not in data:
        log.warning("No Reactome results for '%s'", pathway_query)
        return []

    # Get first pathway hit
    pathway_id = None
    for group in data.get("results", []):
        for entry in group.get("entries", []):
            pathway_id = entry.get("stId")
            pathway_name = entry.get("name", "")
            log.info("Found pathway: %s (%s)", pathway_name, pathway_id)
            break
        if pathway_id:
            break

    if not pathway_id:
        return []

    # Get genes in pathway
    url = f"https://reactome.org/ContentService/data/participants/{pathway_id}"
    participants = http_get_json(url)

    genes = set()
    if participants and isinstance(participants, list):
        for p in participants:
            ref_entities = p.get("refEntities", [])
            for ref in ref_entities:
                gene = ref.get("geneName")
                if gene:
                    genes.update(gene if isinstance(gene, list) else [gene])

    gene_list = sorted(genes)
    log.info("Pathway '%s': %d genes", pathway_query, len(gene_list))
    return gene_list


# ── Pipeline execution ───────────────────────────────────────────


def run_pipeline_for_gene(gene: str, outdir: Path, nextflow_args: list[str] = None) -> bool:
    """Run the Nextflow pipeline for a single gene."""
    gene_dir = outdir / gene
    gene_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "nextflow",
        "run",
        "main.nf",
        "--id",
        gene,
        "--outdir",
        str(gene_dir),
        "-resume",  # resume from cache if available
    ]
    if nextflow_args:
        cmd.extend(nextflow_args)

    log.info("Running pipeline for %s → %s", gene, gene_dir)

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,  # 1 hour max per protein
        )
        if result.returncode == 0:
            log.info("✓ %s completed successfully.", gene)
            return True
        else:
            log.warning(
                "✗ %s failed (exit %d): %s",
                gene,
                result.returncode,
                result.stderr[-500:] if result.stderr else "",
            )
            return False
    except subprocess.TimeoutExpired:
        log.warning("✗ %s timed out after 1 hour.", gene)
        return False
    except FileNotFoundError:
        log.error("Nextflow not found in PATH.")
        return False


# ── Load results ─────────────────────────────────────────────────


def load_protein_results(results_dir: Path) -> dict[str, list[dict]]:
    """
    Load all *_residues.tsv files from a results directory.
    Returns {gene_name: [row_dicts]}
    """
    proteins = {}

    for tsv_path in sorted(results_dir.rglob("*_residues.tsv")):
        gene = tsv_path.stem.replace("_residues", "")
        with open(tsv_path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            proteins[gene] = list(reader)
        log.info("Loaded %s: %d residues", gene, len(proteins[gene]))

    return proteins


# ── Cross-protein comparison ─────────────────────────────────────


def generate_batch_comparison(proteins: dict[str, list[dict]], outdir: Path):
    """
    Generate cross-protein comparative analysis.
    """
    if not proteins:
        log.warning("No proteins to compare.")
        return

    # ── 1. Per-protein summary ───────────────────────────────
    summary_rows = []
    for gene, rows in proteins.items():
        n = len(rows)

        def _avg(field):
            vals = [float(r[field]) for r in rows if r.get(field) and r[field] not in ("", "None")]
            return round(sum(vals) / len(vals), 4) if vals else None

        def _count(field, pred=None):
            if pred:
                return sum(1 for r in rows if pred(r))
            return sum(
                1 for r in rows if r.get(field) and r[field] not in ("", "False", "0", "None")
            )

        # Designability index stats
        di_vals = [
            float(r.get("designability_index", 0)) for r in rows if r.get("designability_index")
        ]
        di_avg = round(sum(di_vals) / len(di_vals), 4) if di_vals else None
        di_engineer = sum(1 for v in di_vals if v >= 0.7)
        di_avoid = sum(1 for v in di_vals if v < 0.15)

        summary_rows.append(
            {
                "gene": gene,
                "length": n,
                "avg_conservation": _avg("conservation_score"),
                "avg_designability_index": di_avg,
                "n_engineer_positions": di_engineer,
                "n_avoid_positions": di_avoid,
                "pct_engineerable": round(di_engineer / n * 100, 1) if n else 0,
                "n_ptm": _count("ptm"),
                "n_disease_variants": _count("disease_association"),
                "n_cosmic_positions": _count(
                    "cosmic_mutation_count",
                    lambda r: (int(r.get("cosmic_mutation_count") or 0)) > 0,
                ),
                "n_cancer_hotspots": _count(
                    "is_cancer_hotspot",
                    lambda r: r.get("is_cancer_hotspot") in ("True", "true", "1"),
                ),
                "n_epitope_residues": _count(
                    "is_epitope", lambda r: r.get("is_epitope") in ("True", "true", "1")
                ),
                "n_pocket_residues": _count(
                    "in_pocket", lambda r: r.get("in_pocket") in ("True", "true", "1")
                ),
                "n_interface_residues": _count(
                    "in_interface", lambda r: r.get("in_interface") in ("True", "true", "1")
                ),
                "n_functional_sites": _count("functional_site"),
                "n_domains": len(set(r.get("domain", "") for r in rows if r.get("domain"))),
            }
        )

    # Write comparison
    comp_path = outdir / "batch_comparison.tsv"
    with open(comp_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)
    log.info("Wrote batch comparison: %s (%d proteins)", comp_path, len(summary_rows))

    # ── 2. Engineering candidates (top DI across all proteins) ───
    candidates = []
    for gene, rows in proteins.items():
        for r in rows:
            di = r.get("designability_index")
            if di and di not in ("", "None"):
                try:
                    di_val = float(di)
                except ValueError:
                    continue
                if di_val >= 0.5:  # only include reasonable candidates
                    candidates.append(
                        {
                            "gene": gene,
                            "position": r.get("position"),
                            "amino_acid": r.get("amino_acid"),
                            "designability_index": di_val,
                            "di_class": r.get("di_class", ""),
                            "conservation_score": r.get("conservation_score", ""),
                            "esm_llr": r.get("esm_llr", ""),
                            "accessibility": r.get("accessibility", ""),
                            "domain": r.get("domain", ""),
                            "ptm": r.get("ptm", ""),
                            "disease_association": r.get("disease_association", ""),
                            "di_rationale": r.get("di_rationale", ""),
                            "esm_top_substitutions": r.get("esm_top_substitutions", ""),
                        }
                    )

    candidates.sort(key=lambda x: -x["designability_index"])
    top_candidates = candidates[:500]  # top 500

    cand_path = outdir / "engineering_candidates.tsv"
    if top_candidates:
        with open(cand_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(top_candidates[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(top_candidates)
        log.info("Wrote engineering candidates: %s (%d positions)", cand_path, len(top_candidates))

    # ── 3. Domain conservation comparison ────────────────────
    domain_stats = defaultdict(
        lambda: {"genes": set(), "positions": 0, "avg_cons": [], "avg_di": []}
    )
    for gene, rows in proteins.items():
        for r in rows:
            domain = r.get("domain", "")
            if not domain:
                continue
            # May have multiple domains semicolon-separated
            for d in domain.split(";"):
                d = d.strip()
                if d:
                    domain_stats[d]["genes"].add(gene)
                    domain_stats[d]["positions"] += 1
                    cons = r.get("conservation_score")
                    if cons and cons not in ("", "None"):
                        domain_stats[d]["avg_cons"].append(float(cons))
                    di = r.get("designability_index")
                    if di and di not in ("", "None"):
                        domain_stats[d]["avg_di"].append(float(di))

    if domain_stats:
        domain_rows = []
        for domain, stats in sorted(domain_stats.items(), key=lambda x: -len(x[1]["genes"])):
            if len(stats["genes"]) < 2:
                continue  # only shared domains
            domain_rows.append(
                {
                    "domain": domain,
                    "n_proteins": len(stats["genes"]),
                    "proteins": "; ".join(sorted(stats["genes"])),
                    "total_positions": stats["positions"],
                    "avg_conservation": round(sum(stats["avg_cons"]) / len(stats["avg_cons"]), 4)
                    if stats["avg_cons"]
                    else "",
                    "avg_designability": round(sum(stats["avg_di"]) / len(stats["avg_di"]), 4)
                    if stats["avg_di"]
                    else "",
                }
            )

        if domain_rows:
            dom_path = outdir / "shared_domains.tsv"
            with open(dom_path, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(domain_rows[0].keys()), delimiter="\t")
                w.writeheader()
                w.writerows(domain_rows)
            log.info("Wrote shared domains: %s (%d domains)", dom_path, len(domain_rows))


# =====================================================================
#  ENTRY POINT
# =====================================================================


def main():
    ap = argparse.ArgumentParser(
        description="Batch protein characterization with cross-protein comparison"
    )
    ap.add_argument("--genes", nargs="+", help="Gene symbols to analyze")
    ap.add_argument("--gene-list", help="File with one gene symbol per line")
    ap.add_argument("--pathway", help="Reactome pathway name to fetch genes from")
    ap.add_argument("--results-dir", help="Compare existing results (skip pipeline)")
    ap.add_argument("--outdir", default="batch_results")
    ap.add_argument(
        "--nextflow-args",
        nargs="*",
        default=[],
        help="Additional Nextflow arguments (e.g. -profile docker)",
    )
    ap.add_argument(
        "--skip-pipeline", action="store_true", help="Skip pipeline execution, only do comparison"
    )
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Determine gene list ──────────────────────────────────
    genes = []

    if args.genes:
        genes = args.genes
    elif args.gene_list:
        with open(args.gene_list) as f:
            genes = [line.strip() for line in f if line.strip() and not line.startswith("#")]
    elif args.pathway:
        genes = fetch_pathway_genes(args.pathway)
    elif args.results_dir:
        # Just compare existing results
        pass
    else:
        ap.error("Provide --genes, --gene-list, --pathway, or --results-dir")

    if genes:
        log.info("Batch analysis for %d proteins: %s", len(genes), ", ".join(genes[:10]))

    # ── Run pipelines ────────────────────────────────────────
    if genes and not args.skip_pipeline:
        successful = []
        failed = []
        for gene in genes:
            ok = run_pipeline_for_gene(gene, outdir, args.nextflow_args)
            (successful if ok else failed).append(gene)

        log.info("Pipeline complete: %d/%d successful.", len(successful), len(genes))
        if failed:
            log.warning("Failed: %s", ", ".join(failed))

    # ── Load and compare ─────────────────────────────────────
    results_path = Path(args.results_dir) if args.results_dir else outdir
    proteins = load_protein_results(results_path)

    if len(proteins) < 1:
        log.error("No protein results found in %s", results_path)
        sys.exit(1)

    generate_batch_comparison(proteins, outdir)

    log.info("Batch analysis complete. Results in %s/", outdir)
    log.info("  batch_comparison.tsv       – per-protein summary")
    log.info("  engineering_candidates.tsv – top mutable positions")
    log.info("  shared_domains.tsv         – cross-protein domain analysis")


if __name__ == "__main__":
    main()
