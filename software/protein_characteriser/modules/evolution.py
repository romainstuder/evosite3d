#!/usr/bin/env python3
"""
modules/evolution.py
=====================
Per-residue evolutionary annotations from Ensembl Compara:

  1. Conservation score  – fraction of orthologs with identical AA
     (from Ensembl Compara pairwise homology endpoint)
  2. Last common ancestor – deepest taxonomic group where AA is conserved
  3. Multiple sequence alignment – the phylogeny-aware protein MSA from
     the Ensembl Compara gene tree (NOT reconstructed from pairwise data)

The MSA is fetched via the Ensembl REST genetree/member/id endpoint,
which returns the real Compara protein tree alignment built with their
pipeline (TreeBeST / phylogeny-guided multiple alignment). This is a
genuine MSA where all species are aligned simultaneously, not a
projection of pairwise alignments.

Inputs:  resolved_ids.json, sequence.fasta
Outputs: evolution.tsv
         orthologs.fasta          (Compara gene tree protein MSA)
         jalview_features.txt
         jalview_annotations.txt
"""

import argparse
import csv
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from taxonomy import tax_bucket_rank, tax_classify

from lib.ensembl import fetch_compara_msa, fetch_pairwise_ortholog

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


# =====================================================================
#  PART 1: Conservation scores (from Ensembl pairwise homology)
# =====================================================================


def compute_conservation(ensembl_gene: str, sequence: str) -> list[dict]:
    """Compute per-residue conservation from Ensembl Compara pairwise orthologs.

    For each position in the human protein sequence, calculates the
    fraction of vertebrate orthologs with an identical amino acid and
    identifies the most distant taxonomic group (last common ancestor)
    where the residue is conserved.

    Args:
        ensembl_gene: Ensembl gene stable ID (e.g. ``"ENSG00000141510"``).
        sequence: Human protein sequence (one-letter amino acid codes).

    Returns:
        List of dicts, one per residue, with keys ``position`` (1-based),
        ``conservation_score`` (float 0--1 or ``None``), and
        ``last_common_ancestor`` (taxonomy bucket name or ``None``).
    """
    seq_len = len(sequence)
    results = [
        {
            "position": i + 1,
            "aa": sequence[i],
            "conservation_score": None,
            "last_common_ancestor": None,
        }
        for i in range(seq_len)
    ]

    data = fetch_pairwise_ortholog(ensembl_gene)

    if not data or "data" not in data or not data["data"]:
        log.warning("No homology data for %s", ensembl_gene)
        return results

    homologies = data["data"][0].get("homologies", [])
    log.info("Received %d ortholog pairs.", len(homologies))
    if not homologies:
        return results

    pos_total = defaultdict(int)
    pos_match = defaultdict(int)
    pos_buckets: dict[int, set] = defaultdict(set)

    for hom in homologies:
        src = hom.get("source", {})
        tgt = hom.get("target", {})
        src_align = src.get("align_seq", "")
        tgt_align = tgt.get("align_seq", "")

        if not src_align or not tgt_align or len(src_align) != len(tgt_align):
            continue

        species = tgt.get("species", "")
        taxon_level = hom.get("taxonomy_level", "")
        bucket = tax_classify(species, taxon_level)

        human_pos = 0
        for col in range(len(src_align)):
            if src_align[col] == "-":
                continue
            human_pos += 1
            if tgt_align[col] == "-":
                continue
            pos_total[human_pos] += 1
            if src_align[col] == tgt_align[col]:
                pos_match[human_pos] += 1
                if bucket:
                    pos_buckets[human_pos].add(bucket)

    for r in results:
        p = r["position"]
        total = pos_total.get(p, 0)
        if total > 0:
            r["conservation_score"] = round(pos_match.get(p, 0) / total, 4)
            buckets = pos_buckets.get(p, set())
            if buckets:
                r["last_common_ancestor"] = min(buckets, key=tax_bucket_rank)
            else:
                r["last_common_ancestor"] = "Human-specific"

    scored = sum(1 for r in results if r["conservation_score"] is not None)
    log.info("Conservation computed for %d / %d positions.", scored, seq_len)
    return results


# =====================================================================
#  PART 3: Jalview output files
# =====================================================================


def write_jalview_features(results: list[dict], seq_id: str, outpath: Path):
    """Write a Jalview features file with conservation tiers and LCA annotations.

    Produces colour-coded features for three conservation tiers (high > 0.9,
    medium > 0.5, low <= 0.5) and a separate group for last-common-ancestor
    labels. The output follows the Jalview features file format.

    Args:
        results: Per-residue conservation dicts as returned by
            :func:`compute_conservation`.
        seq_id: Sequence identifier in the MSA to attach features to.
        outpath: Destination file path.
    """
    with open(outpath, "w") as f:
        f.write("conservation_high\t#d73027\n")
        f.write("conservation_med\t#fee08b\n")
        f.write("conservation_low\t#1a9850\n")

        f.write("\nSTARTGROUP\tConservation\n")
        for r in results:
            pos = r["position"]
            cons = r.get("conservation_score")
            if cons is None:
                continue
            if cons > 0.9:
                ftype, desc = "conservation_high", f"Highly conserved ({cons:.3f})"
            elif cons > 0.5:
                ftype, desc = "conservation_med", f"Moderately conserved ({cons:.3f})"
            else:
                ftype, desc = "conservation_low", f"Low conservation ({cons:.3f})"
            f.write(f"{desc}\t{seq_id}\t-1\t{pos}\t{pos}\t{ftype}\n")
        f.write("ENDGROUP\tConservation\n")

        f.write("\nSTARTGROUP\tLast Common Ancestor\n")
        for r in results:
            pos = r["position"]
            lca = r.get("last_common_ancestor")
            if lca:
                f.write(f"LCA: {lca}\t{seq_id}\t-1\t{pos}\t{pos}\tLCA\n")
        f.write("ENDGROUP\tLast Common Ancestor\n")

    log.info("Wrote Jalview features: %s", outpath)


def write_jalview_annotations(results: list[dict], seq_id: str, outpath: Path):
    """Write a Jalview annotation file with a conservation bar graph.

    Each residue is represented as a coloured bar (red for highly
    conserved, green for variable) in the Jalview ``BAR_GRAPH`` format.

    Args:
        results: Per-residue conservation dicts as returned by
            :func:`compute_conservation`.
        seq_id: Sequence identifier in the MSA to attach annotations to.
        outpath: Destination file path.
    """
    with open(outpath, "w") as f:
        f.write("JALVIEW_ANNOTATION\n")
        f.write("# Conservation from Ensembl Compara pairwise orthologs\n")
        f.write("# MSA from Ensembl Compara gene tree (phylogeny-aware)\n\n")
        f.write(f"SEQUENCE_REF\t{seq_id}\n")

        vals = []
        for r in results:
            cons = r.get("conservation_score")
            if cons is not None:
                if cons > 0.9:
                    color = "ff0000"
                elif cons > 0.7:
                    color = "ff8800"
                elif cons > 0.5:
                    color = "ffcc00"
                elif cons > 0.3:
                    color = "88cc00"
                else:
                    color = "00aa00"
                vals.append(f"{cons:.3f},{cons:.2f},{color}")
            else:
                vals.append("")

        f.write("BAR_GRAPH\tConservation\tPairwise ortholog conservation (0-1)\t")
        f.write("|".join(vals))
        f.write("\n")

    log.info("Wrote Jalview annotations: %s", outpath)


# =====================================================================
#  ENTRY POINT
# =====================================================================


def main():
    """Run the evolutionary analysis pipeline.

    Parses command-line arguments, computes per-residue conservation
    scores from Ensembl Compara pairwise orthologs, optionally fetches
    the Compara gene-tree MSA, and writes output files (TSV, FASTA,
    Jalview features and annotations).
    """
    ap = argparse.ArgumentParser(
        description="Evolutionary analysis: conservation + Compara gene tree MSA"
    )
    ap.add_argument("--ids-json", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument(
        "--skip-msa", action="store_true", help="Skip fetching the Compara gene tree MSA"
    )
    args = ap.parse_args()

    with open(args.ids_json) as f:
        ids = json.load(f)

    record = SeqIO.read(args.fasta, "fasta")
    sequence = str(record.seq)
    outdir = Path(args.outdir)

    ens_prot = ids.get("ensembl_protein")
    ens_gene = ids.get("ensembl_gene")
    if not ens_gene:
        log.error("No Ensembl gene ID – cannot compute conservation.")
        sys.exit(1)

    # ── Part 1: Conservation from pairwise homology ──────────────
    rows = compute_conservation(ens_gene, sequence)

    outpath = outdir / "evolution.tsv"
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["position", "aa", "conservation_score", "last_common_ancestor"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows)
    log.info("Wrote %s", outpath)

    # ── Part 2: Fetch Compara gene tree MSA ──────────────────────
    if args.skip_msa:
        log.info("--skip-msa: skipping gene tree MSA.")
        return

    msa_path = outdir / "orthologs.fasta"
    human_id = fetch_compara_msa(ens_gene, ens_prot, msa_path)

    # ── Part 3: Jalview files ────────────────────────────────────
    if human_id and msa_path.exists():
        write_jalview_features(rows, human_id, outdir / "jalview_features.txt")
        write_jalview_annotations(rows, human_id, outdir / "jalview_annotations.txt")
    elif msa_path.exists():
        # MSA exists but couldn't identify human sequence – still useful
        log.warning(
            "Could not identify human sequence in MSA. Jalview features will use a generic ID."
        )
        fallback_id = f"Homo_sapiens__{ens_prot}"
        write_jalview_features(rows, fallback_id, outdir / "jalview_features.txt")
        write_jalview_annotations(rows, fallback_id, outdir / "jalview_annotations.txt")
    else:
        log.warning("No MSA produced – Jalview files not generated.")


if __name__ == "__main__":
    main()
