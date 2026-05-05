#!/usr/bin/env python3
"""
modules/structure.py
=====================
Per-residue structural annotations:

  1. PDB mapping    – best structure(s) covering each UniProt residue (SIFTS)
  2. Cα coordinates – x, y, z extracted with BioPython's PDB parser
  3. Pocket         – ligand-binding pocket residues from PDBe-KB
  4. Interface      – protein–protein interaction interface from PDBe-KB

Uses BioPython PDBParser / MMCIFParser for coordinate extraction, which
is much more reliable than REST-based atom retrieval.

Inputs:  resolved_ids.json, sequence.fasta
Outputs: structure.tsv
"""

import argparse
import csv
import json
import logging
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.PDB import MMCIFParser, PDBList
from Bio.PDB.Polypeptide import is_aa

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from api_utils import PDBE_GRAPH, PDBE_REST, http_get_json

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


# ── 1. SIFTS: UniProt ↔ PDB residue mapping ─────────────────────


def get_sifts_mappings(uniprot: str) -> list[dict]:
    url = f"{PDBE_REST}/api/mappings/uniprot/{uniprot}"
    data = http_get_json(url)
    if not data or uniprot not in data:
        return []
    records = []
    for pdb_id, info in data[uniprot].get("PDB", {}).items():
        for m in info:
            records.append(
                {
                    "pdb_id": pdb_id,
                    "chain_id": m.get("chain_id"),
                    "struct_asym_id": m.get("struct_asym_id"),
                    "unp_start": m.get("unp_start"),
                    "unp_end": m.get("unp_end"),
                }
            )
    # Sort by coverage width, widest first
    records.sort(key=lambda r: (r["unp_end"] or 0) - (r["unp_start"] or 0), reverse=True)
    return records


def pick_best_structures(mappings: list[dict], max_n: int = 5) -> list[dict]:
    """Greedy set-cover: pick structures that maximize residue coverage."""
    covered = set()
    selected = []
    for m in mappings:
        span = set(range(m["unp_start"] or 0, (m["unp_end"] or 0) + 1))
        new = span - covered
        if new:
            selected.append(m)
            covered |= span
        if len(selected) >= max_n:
            break
    return selected


def build_residue_map(uniprot: str, pdb_id: str, chain_id: str) -> dict[int, int]:
    """SIFTS segment mapping: {uniprot_pos → pdb_residue_number}."""
    url = f"{PDBE_REST}/api/mappings/uniprot_segments/{pdb_id}"
    data = http_get_json(url)
    if not data or pdb_id.lower() not in data:
        return {}
    mapping = {}
    segments = data[pdb_id.lower()].get("UniProt", {}).get(uniprot, {}).get("mappings", [])
    for seg in segments:
        if seg.get("chain_id") != chain_id:
            continue
        us = seg.get("unp_start", 0)
        ue = seg.get("unp_end", 0)
        ps = seg.get("start", {}).get("residue_number", 0)
        if us and ue and ps:
            for off in range(ue - us + 1):
                mapping[us + off] = ps + off
    return mapping


# ── 2. Cα extraction with BioPython ─────────────────────────────


def fetch_ca_coords(pdb_id: str, chain_id: str, tmpdir: str) -> dict[int, tuple]:
    """
    Download structure via BioPython PDBList, parse with MMCIFParser,
    return {residue_number: (x, y, z)} for Cα atoms in the target chain.
    """
    coords = {}
    pdbl = PDBList(verbose=False)

    try:
        # Download mmCIF (preferred modern format)
        path = pdbl.retrieve_pdb_file(pdb_id, pdir=tmpdir, file_format="mmCif")
        if not path or not Path(path).exists():
            log.warning("Could not download %s", pdb_id)
            return coords

        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(pdb_id, path)
    except Exception as e:
        log.warning("BioPython parse error for %s: %s", pdb_id, e)
        return coords

    # Walk the structure: model 0 → target chain → residues → CA atom
    model = structure[0]
    for chain in model:
        if chain.id != chain_id:
            continue
        for residue in chain:
            if not is_aa(residue, standard=True):
                continue
            resnum = residue.id[1]
            if "CA" in residue:
                ca = residue["CA"]
                v = ca.get_vector()
                coords[resnum] = (round(v[0], 3), round(v[1], 3), round(v[2], 3))

    log.info("Extracted %d Cα coords from %s chain %s", len(coords), pdb_id, chain_id)
    return coords


# ── 3. Pocket residues (PDBe-KB) ────────────────────────────────


def fetch_pocket_residues(uniprot: str) -> set[int]:
    positions = set()
    url = f"{PDBE_GRAPH}/uniprot/annotations/{uniprot}"
    data = http_get_json(url)
    if not data or uniprot not in data:
        return positions

    for ann in data[uniprot]:
        atype = ann.get("type", "").lower()
        if not any(kw in atype for kw in ("pocket", "binding", "ligand", "active_site", "site")):
            continue
        for region in ann.get("data", []):
            s = region.get("startIndex") or region.get("start")
            e = region.get("endIndex") or region.get("end")
            if s and e:
                positions.update(range(int(s), int(e) + 1))
            for res in region.get("residues", []):
                p = res.get("startIndex") or res.get("residue_number")
                if p:
                    positions.add(int(p))

    log.info("Pocket residues (PDBe-KB): %d", len(positions))
    return positions


# ── 4. Interface residues (PDBe-KB) ─────────────────────────────


def fetch_interface_residues(uniprot: str) -> set[int]:
    positions = set()
    url = f"{PDBE_GRAPH}/uniprot/interface_residues/{uniprot}"
    data = http_get_json(url)
    if not data or uniprot not in data:
        return positions

    for entry in data[uniprot]:
        if not isinstance(entry, dict):
            continue
        for partner in entry.get("data", []):
            for res in partner.get("residues", []):
                s = res.get("startIndex") or res.get("start")
                e = res.get("endIndex") or res.get("end")
                if s and e:
                    positions.update(range(int(s), int(e) + 1))
                elif s:
                    positions.add(int(s))

    log.info("Interface residues (PDBe-KB): %d", len(positions))
    return positions


# ── Pipeline ─────────────────────────────────────────────────────


def annotate_structure(uniprot: str, sequence: str) -> list[dict]:
    seq_len = len(sequence)
    results = [
        {
            "position": i + 1,
            "aa": sequence[i],
            "pdb_id": None,
            "x": None,
            "y": None,
            "z": None,
            "in_pocket": False,
            "in_interface": False,
        }
        for i in range(seq_len)
    ]

    if not uniprot:
        log.warning("No UniProt – skipping structure.")
        return results

    # Pocket + interface (UniProt-level, independent of PDB choice)
    pockets = fetch_pocket_residues(uniprot)
    interfaces = fetch_interface_residues(uniprot)
    for r in results:
        r["in_pocket"] = r["position"] in pockets
        r["in_interface"] = r["position"] in interfaces

    # SIFTS PDB mapping
    log.info("Fetching SIFTS mappings for %s …", uniprot)
    sifts = get_sifts_mappings(uniprot)
    if not sifts:
        log.warning("No PDB structures for %s", uniprot)
        return results

    best = pick_best_structures(sifts)
    log.info(
        "Selected %d structure(s): %s",
        len(best),
        ", ".join(f"{b['pdb_id']}_{b['chain_id']}" for b in best),
    )

    # Build combined uniprot_pos → (pdb_id, chain, pdb_res)
    unp_to_struct: dict[int, tuple] = {}
    for s in best:
        rmap = build_residue_map(uniprot, s["pdb_id"], s["chain_id"])
        for upos, pres in rmap.items():
            if upos not in unp_to_struct:
                unp_to_struct[upos] = (s["pdb_id"], s["chain_id"], pres)

    # Fetch coordinates with BioPython
    pdb_chains = {(pid, cid) for pid, cid, _ in unp_to_struct.values()}
    all_coords: dict[tuple, dict] = {}

    with tempfile.TemporaryDirectory(prefix="pdb_") as tmpdir:
        for pdb_id, chain_id in pdb_chains:
            log.info("Downloading & parsing %s chain %s …", pdb_id, chain_id)
            all_coords[(pdb_id, chain_id)] = fetch_ca_coords(pdb_id, chain_id, tmpdir)

    # Assign
    for r in results:
        pos = r["position"]
        if pos in unp_to_struct:
            pdb_id, chain_id, pdb_res = unp_to_struct[pos]
            r["pdb_id"] = f"{pdb_id}_{chain_id}"
            xyz = all_coords.get((pdb_id, chain_id), {}).get(pdb_res)
            if xyz:
                r["x"], r["y"], r["z"] = xyz

    covered = sum(1 for r in results if r["pdb_id"])
    with_xyz = sum(1 for r in results if r["x"] is not None)
    log.info("Coverage: %d/%d mapped, %d with Cα coords.", covered, seq_len, with_xyz)
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
    sequence = str(record.seq)

    rows = annotate_structure(ids.get("uniprot"), sequence)

    outpath = Path(args.outdir) / "structure.tsv"
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["position", "aa", "pdb_id", "x", "y", "z", "in_pocket", "in_interface"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows)
    log.info("Wrote %s", outpath)


if __name__ == "__main__":
    main()
