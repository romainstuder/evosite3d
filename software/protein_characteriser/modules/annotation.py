#!/usr/bin/env python3
"""
modules/annotation.py
======================
Per-residue structural & functional context annotations:

  ┌──────────────────────────────────────────────────────────────┐
  │  1. SECONDARY STRUCTURE (DSSP via BioPython)                 │
  │     Helix (H), Sheet (E), Coil (C) assignment from PDB      │
  │     structures using Bio.PDB.DSSP.                           │
  │                                                              │
  │  2. SOLVENT ACCESSIBILITY                                    │
  │     Relative solvent accessibility (RSA) from DSSP.          │
  │     Classified: buried (<0.2), intermediate (0.2-0.4),       │
  │     exposed (>0.4).                                          │
  │                                                              │
  │  3. PROTEIN DOMAINS (InterPro/Pfam)                          │
  │     Domain and family assignments mapped to each residue     │
  │     from InterPro via the Ensembl REST or EBI InterPro API.  │
  │                                                              │
  │  4. ALPHAFOLD pLDDT                                          │
  │     Per-residue confidence score from AlphaFold DB.          │
  │     >90 very high, 70-90 confident, 50-70 low, <50 disorder │
  │                                                              │
  │  5. POPULATION VARIATION (gnomAD)                            │
  │     Per-position allele frequencies from gnomAD via Ensembl, │
  │     distinguishing common (>1%) vs rare variation.           │
  │                                                              │
  │  6. UNIPROT REGIONS                                          │
  │     Signal peptide, transmembrane, coiled-coil, disordered   │
  │     regions from UniProt features.                           │
  └──────────────────────────────────────────────────────────────┘

Inputs:  resolved_ids.json, sequence.fasta
Outputs: annotation.tsv
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
from Bio.PDB import MMCIFParser, PDBList
from Bio.PDB.Polypeptide import is_aa

# DSSP import – try the modern location first
try:
    from Bio.PDB.DSSP import DSSP

    HAS_DSSP = True
except ImportError:
    HAS_DSSP = False

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from api_utils import PDBE_REST, UNIPROT_REST, http_get_json, http_get_text

from lib.ensembl import fetch_translation_overlap

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


# ── 1 & 2. Secondary structure + solvent accessibility (DSSP) ────


def _get_best_pdb(uniprot: str) -> tuple[str, str] | None:
    """Get the best-resolution PDB structure + chain for DSSP."""
    url = f"{PDBE_REST}/api/mappings/best_structures/{uniprot}"
    data = http_get_json(url)
    if not data or uniprot not in data:
        return None
    best = data[uniprot]
    if not best:
        return None
    entry = best[0]  # highest coverage & resolution
    return entry.get("pdb_id"), entry.get("chain_id")


def _get_sifts_residue_map(uniprot: str, pdb_id: str, chain_id: str) -> dict[int, int]:
    """UniProt pos → PDB residue number."""
    url = f"{PDBE_REST}/api/mappings/uniprot_segments/{pdb_id}"
    data = http_get_json(url)
    if not data or pdb_id.lower() not in data:
        return {}
    mapping = {}
    segments = data[pdb_id.lower()].get("UniProt", {}).get(uniprot, {}).get("mappings", [])
    for seg in segments:
        if seg.get("chain_id") != chain_id:
            continue
        us, ue = seg.get("unp_start", 0), seg.get("unp_end", 0)
        ps = seg.get("start", {}).get("residue_number", 0)
        if us and ue and ps:
            for off in range(ue - us + 1):
                mapping[us + off] = ps + off
    return mapping


def compute_dssp(uniprot: str) -> dict[int, dict]:
    """Run DSSP on the best PDB structure, return
    {uniprot_pos: {"secondary_structure": "H"/"E"/"C",
                   "rsa": float,
                   "accessibility": "buried"/"intermediate"/"exposed"}}
    """
    dssp_results: dict[int, dict] = {}

    if not HAS_DSSP:
        log.warning(
            "Bio.PDB.DSSP not available – skipping DSSP analysis. "
            "Install mkdssp (conda install -c salilab dssp)."
        )
        return dssp_results

    best = _get_best_pdb(uniprot)
    if not best:
        log.warning("No PDB structure for DSSP for %s", uniprot)
        return dssp_results

    pdb_id, chain_id = best
    residue_map = _get_sifts_residue_map(uniprot, pdb_id, chain_id)
    if not residue_map:
        log.warning("No SIFTS mapping for %s %s_%s", uniprot, pdb_id, chain_id)
        return dssp_results

    # Reverse map: pdb_res → uniprot_pos
    pdb_to_unp = {v: k for k, v in residue_map.items()}

    with tempfile.TemporaryDirectory(prefix="dssp_") as tmpdir:
        pdbl = PDBList(verbose=False)
        try:
            path = pdbl.retrieve_pdb_file(pdb_id, pdir=tmpdir, file_format="mmCif")
            if not path or not Path(path).exists():
                return dssp_results

            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure(pdb_id, path)
            model = structure[0]

            dssp = DSSP(model, path, dssp="mkdssp")
        except Exception as e:
            log.warning("DSSP failed for %s: %s", pdb_id, e)
            return dssp_results

        # DSSP results: keyed by (chain_id, (' ', resnum, ' '))
        # Use property_dict directly to avoid _translate_id unpacking issues
        for dssp_key, res_info in dssp.property_dict.items():
            # chain = dssp_key[0]  # type: ignore
            # res_info: (index, AA, sec_struct, RSA, ...)
            if len(res_info) < 4:
                continue

            # aa = res_info[1]
            ss = res_info[2]  # H, B, E, G, I, T, S, -
            rsa = res_info[3]

            # Simplify SS: H/G/I → H, E/B → E, rest → C
            if ss in ("H", "G", "I"):
                ss_simple = "H"
            elif ss in ("E", "B"):
                ss_simple = "E"
            else:
                ss_simple = "C"

            # Classify accessibility
            if rsa < 0.2:
                acc_class = "buried"
            elif rsa < 0.4:
                acc_class = "intermediate"
            else:
                acc_class = "exposed"

            # Map back to UniProt position
            # dssp_key format varies; extract residue number
            try:
                if isinstance(dssp_key, tuple) and len(dssp_key) >= 2:
                    res_id = dssp_key[1]
                    if isinstance(res_id, tuple):
                        resnum = res_id[1]
                    else:
                        resnum = int(res_id)
                else:
                    continue
            except (ValueError, TypeError, IndexError):
                continue

            unp_pos = pdb_to_unp.get(resnum)
            if unp_pos:
                dssp_results[unp_pos] = {
                    "secondary_structure": ss_simple,
                    "rsa": round(float(rsa), 4),
                    "accessibility": acc_class,
                }

    log.info("DSSP: %d positions annotated from %s_%s", len(dssp_results), pdb_id, chain_id)
    return dssp_results


def _compute_dssp_alphafold(uniprot: str) -> dict[int, dict]:
    """Run DSSP on the AlphaFold predicted structure for full-length coverage.

    Returns the same format as :func:`compute_dssp`:
    ``{uniprot_pos: {"secondary_structure", "rsa", "accessibility"}}``.
    AlphaFold residue numbering matches UniProt positions directly.
    """
    dssp_results: dict[int, dict] = {}

    if not HAS_DSSP:
        return dssp_results

    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"
    data = http_get_json(url)
    if not data or not isinstance(data, list) or len(data) == 0:
        return dssp_results

    cif_url = data[0].get("cifUrl") or data[0].get("pdbUrl")
    if not cif_url:
        return dssp_results

    with tempfile.TemporaryDirectory(prefix="af_dssp_") as tmpdir:
        cif_text = http_get_text(cif_url)
        if not cif_text:
            return dssp_results

        cif_path = Path(tmpdir) / "af_model.cif"
        cif_path.write_text(cif_text)

        try:
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("AF", str(cif_path))
            model = structure[0]
            dssp = DSSP(model, str(cif_path), dssp="mkdssp")
        except Exception as e:
            log.warning("DSSP on AlphaFold structure failed: %s", e)
            return dssp_results

        for dssp_key, res_info in dssp.property_dict.items():
            if len(res_info) < 4:
                continue
            ss = res_info[2]
            rsa = res_info[3]

            if ss in ("H", "G", "I"):
                ss_simple = "H"
            elif ss in ("E", "B"):
                ss_simple = "E"
            else:
                ss_simple = "C"

            if rsa < 0.2:
                acc_class = "buried"
            elif rsa < 0.4:
                acc_class = "intermediate"
            else:
                acc_class = "exposed"

            try:
                res_id = dssp_key[1]
                resnum = res_id[1] if isinstance(res_id, tuple) else int(res_id)
            except (ValueError, TypeError, IndexError):
                continue

            dssp_results[resnum] = {
                "secondary_structure": ss_simple,
                "rsa": round(float(rsa), 4),
                "accessibility": acc_class,
            }

    log.info("DSSP (AlphaFold): %d positions for %s", len(dssp_results), uniprot)
    return dssp_results


# ── 3. Protein domains (InterPro / Pfam) ────────────────────────


def fetch_domains(uniprot: str) -> dict[int, list[str]]:
    """
    Fetch InterPro/Pfam domain annotations per residue.
    Returns {position: ["PF00870:P53 (1-290)", ...]}
    """
    domains: dict[int, list[str]] = defaultdict(list)

    # InterPro REST
    url = f"https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{uniprot}?format=json"
    data = http_get_json(url)

    if data and "results" in data:
        for result in data["results"]:
            entry_acc = result.get("metadata", {}).get("accession", "")
            entry_name = result.get("metadata", {}).get("name", "")
            # entry_type = result.get("metadata", {}).get("type", "")
            db = result.get("metadata", {}).get("source_database", "")

            label = f"{entry_acc}"
            if entry_name:
                label += f":{entry_name}"
            if db:
                label += f" ({db})"

            proteins = result.get("proteins", [])
            for prot in proteins:
                locations = prot.get("entry_protein_locations", [])
                for loc in locations:
                    fragments = loc.get("fragments", [])
                    for frag in fragments:
                        start = frag.get("start")
                        end = frag.get("end")
                        if start and end:
                            for p in range(int(start), int(end) + 1):
                                domains[p].append(label)

    # Fallback: UniProt features
    if not domains:
        url = f"{UNIPROT_REST}/uniprotkb/{uniprot}.json"
        up_data = http_get_json(url)
        if up_data:
            for feat in up_data.get("features", []):
                ftype = feat.get("type", "")
                if ftype not in ("Domain", "Region"):
                    continue
                desc = feat.get("description", ftype)
                loc = feat.get("location", {})
                start = loc.get("start", {}).get("value")
                end = loc.get("end", {}).get("value")
                if start and end:
                    for p in range(int(start), int(end) + 1):
                        domains[p].append(desc)

    log.info("Domain annotations: %d positions", len(domains))
    return domains


# ── 4. AlphaFold pLDDT ──────────────────────────────────────────


def fetch_alphafold_plddt(uniprot: str) -> dict[int, float]:
    """
    Fetch per-residue pLDDT confidence from AlphaFold DB.
    Uses the PAE/pLDDT JSON summary endpoint.
    """
    plddt: dict[int, float] = {}

    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"
    data = http_get_json(url)

    if not data or not isinstance(data, list) or len(data) == 0:
        log.info("No AlphaFold prediction for %s", uniprot)
        return plddt

    entry = data[0]
    plddt_url = entry.get("pdbUrl")  # mmCIF with pLDDT in B-factor column
    cif_url = entry.get("cifUrl") or plddt_url

    if not cif_url:
        return plddt

    # Parse mmCIF B-factors (pLDDT) using BioPython
    with tempfile.TemporaryDirectory(prefix="af_") as tmpdir:
        cif_text = http_get_text(cif_url)
        if not cif_text:
            return plddt

        cif_path = Path(tmpdir) / "af_model.cif"
        cif_path.write_text(cif_text)

        try:
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("AF", str(cif_path))
            model = structure[0]
            for chain in model:
                for residue in chain:
                    if not is_aa(residue, standard=True):
                        continue
                    resnum = residue.id[1]
                    # pLDDT is stored in B-factor of CA atom
                    if "CA" in residue:
                        bfactor = residue["CA"].get_bfactor()
                        plddt[resnum] = round(bfactor, 2)
        except Exception as e:
            log.warning("AlphaFold CIF parse error: %s", e)

    log.info("AlphaFold pLDDT: %d positions for %s", len(plddt), uniprot)
    return plddt


def _classify_plddt(score: float | None) -> str:
    """Classify AlphaFold confidence."""
    if score is None:
        return ""
    if score >= 90:
        return "very_high"
    elif score >= 70:
        return "confident"
    elif score >= 50:
        return "low"
    else:
        return "disordered"


# ── 5. Population variation (gnomAD via Ensembl) ────────────────


def fetch_gnomad_frequencies(ensembl_protein: str) -> dict[int, dict]:
    """
    Fetch per-position gnomAD allele frequencies for missense variants.
    Returns {position: {"max_af": float, "n_variants": int,
                         "frequency_class": "common"/"rare"/"singleton"}}
    """
    gnomad: dict[int, dict] = defaultdict(lambda: {"max_af": 0.0, "n_variants": 0})

    if not ensembl_protein:
        return gnomad

    data = fetch_translation_overlap(ensembl_protein, feature="transcript_variation")
    if not data:
        return gnomad

    for var in data:
        if "missense_variant" not in var.get("consequence_terms", []):
            continue

        start = var.get("start")
        if start is None:
            continue
        pos = int(start)

        # gnomAD frequencies are in the "frequencies" field or
        # individual population fields
        freqs = var.get("frequencies", {})
        if not freqs:
            # Try minor_allele_freq
            maf = var.get("minor_allele_freq")
            if maf is not None:
                try:
                    af = float(maf)
                    gnomad[pos]["max_af"] = max(gnomad[pos]["max_af"], af)
                    gnomad[pos]["n_variants"] += 1
                except (ValueError, TypeError):
                    pass
        else:
            for pop, af_data in freqs.items():
                if isinstance(af_data, dict):
                    af = af_data.get("af") or af_data.get("allele_freq")
                elif isinstance(af_data, (int, float)):
                    af = af_data
                else:
                    continue
                if af is not None:
                    try:
                        gnomad[pos]["max_af"] = max(gnomad[pos]["max_af"], float(af))
                        gnomad[pos]["n_variants"] += 1
                    except (ValueError, TypeError):
                        pass

    # Classify
    for pos in gnomad:
        af = gnomad[pos]["max_af"]
        if af >= 0.01:
            gnomad[pos]["frequency_class"] = "common"
        elif af >= 0.0001:
            gnomad[pos]["frequency_class"] = "rare"
        elif af > 0:
            gnomad[pos]["frequency_class"] = "singleton"
        else:
            gnomad[pos]["frequency_class"] = ""

    log.info("gnomAD: %d positions with frequency data", len(gnomad))
    return gnomad


# ── 6. UniProt region features ──────────────────────────────────


def fetch_uniprot_regions(uniprot: str) -> dict[int, list[str]]:
    """
    Signal peptide, transmembrane, coiled-coil, disordered, etc.
    """
    regions: dict[int, list[str]] = defaultdict(list)

    url = f"{UNIPROT_REST}/uniprotkb/{uniprot}.json"
    data = http_get_json(url)
    if not data:
        return regions

    _REGION_TYPES = {
        "Signal peptide",
        "Transit peptide",
        "Transmembrane",
        "Intramembrane",
        "Topological domain",
        "Coiled coil",
        "Compositional bias",
        "Zinc finger",
        "DNA binding",
        "Nucleotide binding",
        "Motif",
        "Chain",
        "Peptide",
    }

    for feat in data.get("features", []):
        ftype = feat.get("type", "")
        if ftype not in _REGION_TYPES:
            continue
        desc = feat.get("description", ftype)
        label = f"{ftype}: {desc}" if desc and desc != ftype else ftype

        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None:
            continue
        if end is None:
            end = start

        for p in range(int(start), int(end) + 1):
            regions[p].append(label)

    log.info("UniProt regions: %d positions annotated", len(regions))
    return regions


# =====================================================================
#  PUBLIC PIPELINE
# =====================================================================

COLUMNS = [
    "position",
    "aa",
    "secondary_structure",
    "rsa",
    "accessibility",
    "domain",
    "alphafold_plddt",
    "alphafold_confidence",
    "gnomad_max_af",
    "gnomad_frequency_class",
    "uniprot_region",
]


def annotate_annotation(ids: dict, sequence: str) -> list[dict]:
    """Full structural/functional context annotation."""
    seq_len = len(sequence)
    uniprot = ids.get("uniprot")
    ens_prot = ids.get("ensembl_protein")

    # ── Fetch all data ───────────────────────────────────────
    dssp_data = compute_dssp(uniprot) if uniprot else {}
    # Fill gaps with AlphaFold-predicted secondary structure
    if uniprot and len(dssp_data) < seq_len:
        af_dssp = _compute_dssp_alphafold(uniprot)
        for pos, val in af_dssp.items():
            if pos not in dssp_data:
                dssp_data[pos] = val
        log.info("DSSP after AlphaFold gap-fill: %d / %d positions", len(dssp_data), seq_len)
    domain_data = fetch_domains(uniprot) if uniprot else {}
    plddt_data = fetch_alphafold_plddt(uniprot) if uniprot else {}
    gnomad_data = fetch_gnomad_frequencies(ens_prot)
    region_data = fetch_uniprot_regions(uniprot) if uniprot else {}

    # ── Assemble ─────────────────────────────────────────────
    results = []
    for pos in range(1, seq_len + 1):
        row = {"position": pos, "aa": sequence[pos - 1]}

        # DSSP
        dssp = dssp_data.get(pos, {})
        row["secondary_structure"] = dssp.get("secondary_structure", "")
        row["rsa"] = dssp.get("rsa", "")
        row["accessibility"] = dssp.get("accessibility", "")

        # Domains
        doms = domain_data.get(pos, [])
        row["domain"] = "; ".join(doms) if doms else ""

        # AlphaFold
        plddt_val = plddt_data.get(pos)
        row["alphafold_plddt"] = plddt_val if plddt_val is not None else ""
        row["alphafold_confidence"] = _classify_plddt(plddt_val)

        # gnomAD
        gn = gnomad_data.get(pos, {})
        row["gnomad_max_af"] = gn.get("max_af", "") if gn.get("max_af", 0) > 0 else ""
        row["gnomad_frequency_class"] = gn.get("frequency_class", "")

        # UniProt regions
        regs = region_data.get(pos, [])
        row["uniprot_region"] = "; ".join(regs) if regs else ""

        results.append(row)

    n_ss = sum(1 for r in results if r["secondary_structure"])
    n_dom = sum(1 for r in results if r["domain"])
    n_af = sum(1 for r in results if r["alphafold_plddt"] != "")
    log.info("Annotation: %d SS, %d domain, %d AlphaFold positions", n_ss, n_dom, n_af)

    return results


# ── Entry point ──────────────────────────────────────────────────


def main():
    ap = argparse.ArgumentParser(description="Structural & functional context annotations")
    ap.add_argument("--ids-json", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    with open(args.ids_json) as f:
        ids = json.load(f)
    record = SeqIO.read(args.fasta, "fasta")

    rows = annotate_annotation(ids, str(record.seq))

    outpath = Path(args.outdir) / "annotation.tsv"
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    log.info("Wrote %s", outpath)


if __name__ == "__main__":
    main()
