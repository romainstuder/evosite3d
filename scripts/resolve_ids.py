#!/usr/bin/env python3
"""
bin/resolve_ids.py
===================
Resolves a user identifier (gene symbol / UniProt / Ensembl protein)
to canonical cross-references and fetches the protein sequence.

Inputs:  --id <ID>  --id-type gene|uniprot|ensembl_protein
Outputs: resolved_ids.json   (all cross-references)
         sequence.fasta      (canonical protein sequence)

This is the first step in the Nextflow pipeline; downstream processes
depend on its outputs.
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from urllib.parse import quote

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Allow running from project root or as standalone
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from api_utils import ENSEMBL_REST, UNIPROT_REST, http_get_json

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


# ── Internal helpers ─────────────────────────────────────────────


def _ensembl_gene_from_symbol(symbol: str) -> str | None:
    url = f"{ENSEMBL_REST}/xrefs/symbol/homo_sapiens/{quote(symbol)}?external_db=HGNC"
    data = http_get_json(url)
    if data:
        for e in data:
            if e.get("id", "").startswith("ENSG"):
                return e["id"]
    return None


def _canonical_protein_from_gene(ensembl_gene: str) -> str | None:
    url = f"{ENSEMBL_REST}/lookup/id/{ensembl_gene}?expand=1"
    data = http_get_json(url)
    if not data or "Transcript" not in data:
        return None
    for tx in data["Transcript"]:
        if tx.get("is_canonical") == 1:
            prot = tx.get("Translation", {}).get("id")
            if prot:
                return prot
    for tx in data["Transcript"]:
        prot = tx.get("Translation", {}).get("id")
        if prot:
            return prot
    return None


def _uniprot_from_ensembl(ensembl_protein: str) -> str | None:
    url = f"{ENSEMBL_REST}/xrefs/id/{ensembl_protein}?all_levels=1"
    data = http_get_json(url)
    if not data:
        return None
    for e in data:
        if e.get("dbname") == "Uniprot/SWISSPROT":
            return e["primary_id"]
    for e in data:
        if "Uniprot" in e.get("dbname", ""):
            return e["primary_id"]
    return None


def _ensembl_from_uniprot(uniprot: str) -> tuple[str | None, str | None]:
    url = f"{ENSEMBL_REST}/xrefs/symbol/homo_sapiens/{uniprot}?external_db=Uniprot/SWISSPROT"
    data = http_get_json(url)
    ens_gene = None
    if data:
        for e in data:
            if e.get("id", "").startswith("ENSG"):
                ens_gene = e["id"]
                break
    ens_prot = _canonical_protein_from_gene(ens_gene) if ens_gene else None
    return ens_gene, ens_prot


# ── Main resolver ────────────────────────────────────────────────


def resolve(identifier: str, id_type: str) -> dict:
    r = dict(uniprot=None, ensembl_gene=None, ensembl_protein=None, gene_symbol=None, sequence=None)

    if id_type == "gene":
        r["gene_symbol"] = identifier
        r["ensembl_gene"] = _ensembl_gene_from_symbol(identifier)
        if r["ensembl_gene"]:
            r["ensembl_protein"] = _canonical_protein_from_gene(r["ensembl_gene"])
        if r["ensembl_protein"]:
            r["uniprot"] = _uniprot_from_ensembl(r["ensembl_protein"])

    elif id_type == "uniprot":
        r["uniprot"] = identifier
        r["ensembl_gene"], r["ensembl_protein"] = _ensembl_from_uniprot(identifier)

    elif id_type == "ensembl_protein":
        r["ensembl_protein"] = identifier
        r["uniprot"] = _uniprot_from_ensembl(identifier)
        data = http_get_json(f"{ENSEMBL_REST}/lookup/id/{identifier}")
        if data:
            r["ensembl_gene"] = data.get("Parent")

    # Fetch sequence + gene symbol from UniProt
    if r["uniprot"]:
        entry = http_get_json(f"{UNIPROT_REST}/uniprotkb/{r['uniprot']}.json")
        if entry:
            r["sequence"] = entry.get("sequence", {}).get("value")
            if not r["gene_symbol"]:
                genes = entry.get("genes", [])
                if genes:
                    r["gene_symbol"] = genes[0].get("geneName", {}).get("value")

    # Fallback sequence from Ensembl
    if not r["sequence"] and r["ensembl_protein"]:
        data = http_get_json(f"{ENSEMBL_REST}/sequence/id/{r['ensembl_protein']}?type=protein")
        if data:
            r["sequence"] = data.get("seq")

    if not r["sequence"]:
        raise ValueError(f"Cannot resolve sequence for '{identifier}' (id_type={id_type})")

    return r


# ── Entry point ──────────────────────────────────────────────────


def main():
    ap = argparse.ArgumentParser(description="Resolve protein identifiers.")
    ap.add_argument("--id", required=True)
    ap.add_argument("--id-type", default="gene", choices=["gene", "uniprot", "ensembl_protein"])
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    ids = resolve(args.id, args.id_type)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Write JSON
    with open(outdir / "resolved_ids.json", "w") as f:
        json.dump(ids, f, indent=2)
    log.info(
        "Wrote resolved_ids.json  (gene=%s  uniprot=%s  ensembl_prot=%s  len=%d)",
        ids["gene_symbol"],
        ids["uniprot"],
        ids["ensembl_protein"],
        len(ids["sequence"]),
    )

    # Write FASTA using BioPython
    record = SeqRecord(
        Seq(ids["sequence"]),
        id=ids["uniprot"] or ids["ensembl_protein"] or args.id,
        name=ids["gene_symbol"] or args.id,
        description=f"{ids['gene_symbol']} | {ids['uniprot']} | {ids['ensembl_protein']}",
    )
    SeqIO.write(record, outdir / "sequence.fasta", "fasta")
    log.info("Wrote sequence.fasta")


if __name__ == "__main__":
    main()
