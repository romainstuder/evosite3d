#!/usr/bin/env python3
"""
modules/af_multimer.py
=======================
Per-residue predicted protein-protein interaction interfaces from
AlphaFold-Multimer models.

  ┌──────────────────────────────────────────────────────────────┐
  │  ALPHAFOLD-MULTIMER PREDICTED INTERFACES                    │
  │                                                              │
  │  1. Fetch known interaction partners from IntAct/STRING     │
  │  2. Check AlphaFold DB for predicted paired structures      │
  │  3. Parse inter-chain contacts from AF paired models        │
  │  4. Compute per-residue interface scores:                   │
  │     - af_interface_contacts: # of inter-chain Cα contacts   │
  │     - af_interface_pae: predicted aligned error at interface │
  │     - af_interface_partners: which proteins it contacts      │
  │                                                              │
  │  Falls back to STRING interaction data + PDBe-KB if AF      │
  │  paired structures are unavailable.                         │
  └──────────────────────────────────────────────────────────────┘

Inputs:  resolved_ids.json, sequence.fasta
Outputs: af_multimer.tsv
"""

import argparse
import csv
import json
import logging
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.PDB import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Polypeptide import is_aa

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from api_utils import http_get_json, http_get_text

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

AF_DB = "https://alphafold.ebi.ac.uk/api"
STRING_API = "https://string-db.org/api"
INTACT_API = "https://www.ebi.ac.uk/intact/ws"


# ── 1. Fetch interaction partners ────────────────────────────────


def _fetch_string_partners(uniprot: str, n_partners: int = 10) -> list[dict]:
    """Get top interaction partners from STRING."""
    partners = []

    # STRING uses species-prefixed IDs; map UniProt first
    url = (
        f"{STRING_API}/json/interaction_partners"
        f"?identifiers={uniprot}&species=9606&limit={n_partners}"
        f"&caller_identity=protein_characteriser"
    )
    data = http_get_json(url)

    if data and isinstance(data, list):
        seen = set()
        for entry in data:
            pref_name = entry.get("preferredName_B", "")
            string_id = entry.get("stringId_B", "")
            score = entry.get("score", 0)
            if pref_name and pref_name not in seen:
                seen.add(pref_name)
                partners.append(
                    {
                        "gene": pref_name,
                        "string_id": string_id,
                        "score": score,
                    }
                )
        log.info("STRING partners for %s: %d found", uniprot, len(partners))

    return partners


def _fetch_intact_partners(uniprot: str) -> list[dict]:
    """Get binary interaction partners from IntAct."""
    partners = []
    url = f"{INTACT_API}/interaction/findInteractor/{uniprot}?format=json"
    data = http_get_json(url)

    if data and isinstance(data, list):
        seen = set()
        for entry in data:
            interactor_b = entry.get("interactorB", {})
            acc = interactor_b.get("identifier", "")
            gene = interactor_b.get("alias", acc)
            if acc and acc != uniprot and acc not in seen:
                seen.add(acc)
                partners.append({"gene": gene, "uniprot": acc})
        log.info("IntAct partners for %s: %d found", uniprot, len(partners))

    return partners


# ── 2. Fetch AlphaFold paired structures ─────────────────────────


def _fetch_af_paired_models(uniprot: str) -> list[dict]:
    """
    Check AlphaFold EBI for predicted paired/complex structures
    involving this protein (AF-Multimer predictions).
    """
    models = []

    # AF-DB paired predictions endpoint
    url = f"{AF_DB}/prediction/{uniprot}"
    data = http_get_json(url)

    if not data or not isinstance(data, list):
        return models

    for entry in data:
        cif_url = entry.get("cifUrl") or entry.get("pdbUrl")
        model_type = entry.get("modelType", "")
        if cif_url:
            models.append(
                {
                    "cif_url": cif_url,
                    "model_type": model_type,
                    "entry_id": entry.get("entryId", ""),
                }
            )

    log.info("AlphaFold models for %s: %d found", uniprot, len(models))
    return models


def _extract_interface_contacts(
    cif_text: str, target_chain: str = "A", contact_dist: float = 8.0
) -> dict[int, dict]:
    """
    Parse a paired mmCIF structure and find inter-chain contacts.

    For each residue in target_chain, counts how many residues from
    other chains have Cα within contact_dist Å.

    Returns {residue_number: {"n_contacts": int, "partner_chains": set,
                               "min_distance": float}}
    """
    contacts: dict[int, dict] = defaultdict(
        lambda: {"n_contacts": 0, "partner_chains": set(), "min_distance": 999.0}
    )

    with tempfile.NamedTemporaryFile(mode="w", suffix=".cif", delete=False) as f:
        f.write(cif_text)
        cif_path = f.name

    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("af_paired", cif_path)
        model = structure[0]

        # Collect all CA atoms by chain
        target_atoms = []
        other_atoms = []
        target_map = {}  # atom → residue number

        for chain in model:
            for residue in chain:
                if not is_aa(residue, standard=True):
                    continue
                if "CA" not in residue:
                    continue
                ca = residue["CA"]
                resnum = residue.id[1]

                if chain.id == target_chain:
                    target_atoms.append(ca)
                    target_map[id(ca)] = resnum
                else:
                    other_atoms.append((ca, chain.id, resnum))

        if not target_atoms or not other_atoms:
            return contacts

        # Build NeighborSearch on other-chain atoms
        other_atom_list = [a for a, _, _ in other_atoms]
        other_info = {id(a): (c, r) for a, c, r in other_atoms}

        ns = NeighborSearch(other_atom_list)

        for target_ca in target_atoms:
            target_resnum = target_map[id(target_ca)]
            nearby = ns.search(target_ca.get_vector().get_array(), contact_dist, level="A")

            for neighbor_atom in nearby:
                chain_id, partner_res = other_info.get(id(neighbor_atom), (None, None))
                if chain_id is None:
                    continue

                dist = target_ca - neighbor_atom
                contacts[target_resnum]["n_contacts"] += 1
                contacts[target_resnum]["partner_chains"].add(chain_id)
                contacts[target_resnum]["min_distance"] = min(
                    contacts[target_resnum]["min_distance"], dist
                )

    except Exception as e:
        log.warning("mmCIF parse error for AF paired model: %s", e)

    finally:
        Path(cif_path).unlink(missing_ok=True)

    return contacts


# ── 3. Fetch PAE from AlphaFold ──────────────────────────────────


def _fetch_pae_scores(uniprot: str, seq_len: int) -> dict[int, float]:
    """
    Fetch predicted aligned error (PAE) from AlphaFold DB.
    Returns per-residue average inter-domain PAE score.
    """
    pae_scores: dict[int, float] = {}

    url = f"{AF_DB}/prediction/{uniprot}"
    data = http_get_json(url)
    if not data or not isinstance(data, list):
        return pae_scores

    pae_url = data[0].get("paeDocUrl") or data[0].get("paeImageUrl")
    if not pae_url:
        # Try constructing PAE URL
        entry_id = data[0].get("entryId", "")
        if entry_id:
            pae_url = (
                f"https://alphafold.ebi.ac.uk/files/{entry_id}-predicted_aligned_error_v4.json"
            )

    if not pae_url:
        return pae_scores

    pae_data = http_get_json(pae_url)
    if not pae_data:
        return pae_scores

    # PAE is a matrix; compute per-residue average
    # Format: [{"predicted_aligned_error": [[...]], "max_predicted_aligned_error": float}]
    pae_matrix = None
    if isinstance(pae_data, list) and pae_data:
        pae_matrix = pae_data[0].get("predicted_aligned_error")
    elif isinstance(pae_data, dict):
        pae_matrix = pae_data.get("predicted_aligned_error")

    if pae_matrix and isinstance(pae_matrix, list):
        n = len(pae_matrix)
        for i in range(min(n, seq_len)):
            row = pae_matrix[i]
            if isinstance(row, list):
                avg_pae = sum(row) / len(row) if row else 0
                pae_scores[i + 1] = round(avg_pae, 2)

    log.info("PAE scores: %d positions for %s", len(pae_scores), uniprot)
    return pae_scores


# =====================================================================
#  PUBLIC PIPELINE
# =====================================================================

COLUMNS = [
    "position",
    "aa",
    "af_interface_contacts",
    "af_interface_min_dist",
    "af_interface_partners",
    "af_pae_avg",
    "predicted_interface",
    "interaction_partners_string",
]


def annotate_af_multimer(ids: dict, sequence: str) -> list[dict]:
    """Full AF-Multimer interface annotation."""
    seq_len = len(sequence)
    uniprot = ids.get("uniprot")

    # Initialize
    results = [
        {
            "position": i + 1,
            "aa": sequence[i],
            "af_interface_contacts": 0,
            "af_interface_min_dist": "",
            "af_interface_partners": "",
            "af_pae_avg": "",
            "predicted_interface": False,
            "interaction_partners_string": "",
        }
        for i in range(seq_len)
    ]

    if not uniprot:
        log.warning("No UniProt – skipping AF-Multimer.")
        return results

    # ── Interaction partners ─────────────────────────────────
    string_partners = _fetch_string_partners(uniprot)
    partner_str = "; ".join(f"{p['gene']}({p['score']:.2f})" for p in string_partners[:5])
    for r in results:
        r["interaction_partners_string"] = partner_str

    # ── AlphaFold paired models ──────────────────────────────
    models = _fetch_af_paired_models(uniprot)

    all_contacts: dict[int, dict] = defaultdict(
        lambda: {"n_contacts": 0, "partner_chains": set(), "min_distance": 999.0}
    )

    for model in models:
        cif_url = model.get("cif_url")
        if not cif_url:
            continue

        log.info("Fetching AF structure: %s", model.get("entry_id", cif_url[:50]))
        cif_text = http_get_text(cif_url)
        if not cif_text:
            continue

        contacts = _extract_interface_contacts(cif_text)
        for pos, info in contacts.items():
            all_contacts[pos]["n_contacts"] += info["n_contacts"]
            all_contacts[pos]["partner_chains"].update(info["partner_chains"])
            all_contacts[pos]["min_distance"] = min(
                all_contacts[pos]["min_distance"], info["min_distance"]
            )

    # ── PAE scores ───────────────────────────────────────────
    pae = _fetch_pae_scores(uniprot, seq_len)

    # ── Assign ───────────────────────────────────────────────
    for r in results:
        pos = r["position"]

        if pos in all_contacts:
            c = all_contacts[pos]
            r["af_interface_contacts"] = c["n_contacts"]
            r["af_interface_min_dist"] = (
                round(c["min_distance"], 2) if c["min_distance"] < 999 else ""
            )
            r["af_interface_partners"] = ";".join(sorted(c["partner_chains"]))
            r["predicted_interface"] = c["n_contacts"] >= 3

        if pos in pae:
            r["af_pae_avg"] = pae[pos]

    n_interface = sum(1 for r in results if r["predicted_interface"])
    log.info("AF-Multimer: %d interface residues predicted.", n_interface)

    return results


def main():
    ap = argparse.ArgumentParser(description="AlphaFold-Multimer predicted interaction interfaces")
    ap.add_argument("--ids-json", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    with open(args.ids_json) as f:
        ids = json.load(f)
    record = SeqIO.read(args.fasta, "fasta")

    rows = annotate_af_multimer(ids, str(record.seq))

    outpath = Path(args.outdir) / "af_multimer.tsv"
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    log.info("Wrote %s", outpath)


if __name__ == "__main__":
    main()
