#!/usr/bin/env python3
"""
bin/merge_annotations.py
=========================
Joins all 8 per-module TSV outputs (keyed on 'position') with the
protein sequence to produce the final combined residue-level TSV.

Inputs:  sequence.fasta + 8 module TSVs
Outputs: <gene>_residues.tsv
"""

import argparse
import csv
import json
import logging
from pathlib import Path

from Bio import SeqIO

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

COLUMNS = [
    "position",
    "amino_acid",
    # ── Evolution ────────────────────────────────
    "conservation_score",
    "last_common_ancestor",
    # ── Structure (3D) ──────────────────────────
    "pdb_id",
    "x",
    "y",
    "z",
    "in_pocket",
    "in_interface",
    # ── Structural & functional context ──────────
    "secondary_structure",
    "rsa",
    "accessibility",
    "domain",
    "alphafold_plddt",
    "alphafold_confidence",
    "gnomad_max_af",
    "gnomad_frequency_class",
    "uniprot_region",
    # ── Biophysics ──────────────────────────────
    "hydrophobicity",
    "hydrophobicity_window9",
    "charge_ph7",
    "residue_mw",
    "residue_volume",
    "aa_class",
    "disorder_probability",
    "is_disordered",
    "functional_site",
    "elm_motif",
    "n_isoforms",
    "is_alternatively_spliced",
    "codon",
    "is_cpg_codon",
    "is_rare_codon",
    "proteomics_evidence",
    # ── PTM ──────────────────────────────────────
    "ptm",
    # ── Disease (germline) ───────────────────────
    "disease_association",
    "gene_diseases",
    # ── Cancer (somatic) ─────────────────────────
    "cosmic_ids",
    "cosmic_mutation_count",
    "cosmic_aa_changes",
    "cosmic_cancer_types",
    "is_cancer_hotspot",
    "hotspot_score",
    "cancer_tier",
    "mutagenesis_data",
    "gene_cancer_role",
    # ── Designability ────────────────────────────
    "is_epitope",
    "epitope_type",
    "epitope_source",
    "epitope_score",
    "sift_min_score",
    "sift_median_score",
    "polyphen_max_score",
    "mutational_risk",
    "designability_note",
    # ── ESM-2 language model ────────────────────
    "esm_llr",
    "esm_entropy",
    "esm_entropy_norm",
    "esm_wt_prob",
    "esm_top_substitutions",
    "esm_model",
    # ── AF-Multimer predicted interfaces ────────
    "af_interface_contacts",
    "af_interface_min_dist",
    "af_interface_partners",
    "af_pae_avg",
    "predicted_interface",
    "interaction_partners_string",
    # ── Selection pressure ──────────────────────
    "dnds_alpha",
    "dnds_beta",
    "prob_positive_selection",
    "prob_negative_selection",
    "selection_class",
    "human_specific_substitution",
    "human_lineage_change",
    "ancestral_aa",
    "convergent_evolution",
    "convergent_lineages",
    "top_coevo_partners",
    "max_mi_apc",
    "evolutionary_rate",
    # ── Pervasive selection (CodeML) ──────────
    "codeml_beb_prob",
    "codeml_post_mean_omega",
    "codeml_post_se_omega",
    "codeml_estimated_omega",
    "codeml_omega_se",
    "codeml_selection_class",
    "codeml_lrt_m1a_m2a_pval",
    "codeml_lrt_m7_m8_pval",
    "codeml_lrt_m8_m8a_pval",
    "codeml_significant",
]


def load_tsv(path: Path) -> dict[int, dict]:
    """Load a TSV into {position: {col: val, …}}."""
    rows = {}
    if not path.exists():
        log.warning("Missing file: %s – columns will be empty.", path)
        return rows
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pos = int(row.pop("position"))
            rows[pos] = row
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids-json", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--evolution", required=True)
    ap.add_argument("--structure", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--biophysics", required=True)
    ap.add_argument("--ptm", required=True)
    ap.add_argument("--disease", required=True)
    ap.add_argument("--cancer", required=True)
    ap.add_argument("--designability", required=True)
    ap.add_argument("--esm-scoring", required=True)
    ap.add_argument("--af-multimer", required=True)
    ap.add_argument("--selection", required=True)
    ap.add_argument("--pervasive-selection", required=True)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    with open(args.ids_json) as f:
        ids = json.load(f)

    record = SeqIO.read(args.fasta, "fasta")
    sequence = str(record.seq)

    # Load all 10 module outputs
    evo = load_tsv(Path(args.evolution))
    struct = load_tsv(Path(args.structure))
    annot = load_tsv(Path(args.annotation))
    bioph = load_tsv(Path(args.biophysics))
    ptm = load_tsv(Path(args.ptm))
    dis = load_tsv(Path(args.disease))
    cancer = load_tsv(Path(args.cancer))
    design = load_tsv(Path(args.designability))
    esm = load_tsv(Path(args.esm_scoring))
    afm = load_tsv(Path(args.af_multimer))
    sel = load_tsv(Path(args.selection))
    psel = load_tsv(Path(args.pervasive_selection))

    # Merge – order matches COLUMNS grouping
    merged = []
    for i, aa in enumerate(sequence):
        pos = i + 1
        row = {"position": pos, "amino_acid": aa}
        row.update(evo.get(pos, {}))
        row.update(struct.get(pos, {}))
        row.update(annot.get(pos, {}))
        row.update(bioph.get(pos, {}))
        row.update(ptm.get(pos, {}))
        row.update(dis.get(pos, {}))
        row.update(cancer.get(pos, {}))
        row.update(design.get(pos, {}))
        row.update(esm.get(pos, {}))
        row.update(afm.get(pos, {}))
        row.update(sel.get(pos, {}))
        row.update(psel.get(pos, {}))
        merged.append(row)

    gene = ids.get("gene_symbol") or ids.get("uniprot") or "protein"
    outpath = Path(args.outdir) / f"{gene}_residues.tsv"

    with open(outpath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(merged)

    log.info("Wrote %d residues × %d columns → %s", len(merged), len(COLUMNS), outpath)


if __name__ == "__main__":
    main()
