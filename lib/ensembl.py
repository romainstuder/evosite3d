"""
modules/ensembl.py
===================
Shared Ensembl REST API helpers used by evolution, selection, and
pervasive_selection modules.

Provides:
    - ``fetch_cds_sequences`` – CDS nucleotide sequences for protein IDs
    - ``fetch_gene_tree_newick`` – Compara gene tree in Newick format
    - ``fetch_compara_msa`` – phylogeny-aware protein MSA from gene tree
    - ``extract_protein_ids_from_msa`` – parse Ensembl protein IDs from FASTA
"""

import logging
import re
import sys
from pathlib import Path
from urllib.request import Request, urlopen

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
from api_utils import ENSEMBL_REST, http_get_json

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


# =====================================================================
#  MSA helpers
# =====================================================================


def extract_protein_ids_from_msa(msa_path: Path) -> list[str]:
    """Extract Ensembl protein IDs from MSA FASTA headers.

    Parses sequence identifiers looking for patterns like ``ENSP...``,
    ``ENSMUSP...``, etc. Also handles the ``species__accession`` format
    produced by the evolution module.

    Args:
        msa_path: Path to the protein MSA FASTA file.

    Returns:
        List of Ensembl protein stable IDs found in the MSA.
    """
    ids = []
    for record in SeqIO.parse(msa_path, "fasta"):
        match = re.search(r"(ENS\w*P\d+)", record.id)
        if match:
            ids.append(match.group(1))
        else:
            parts = record.id.split("__")
            if len(parts) >= 2:
                candidate = parts[-1]
                if candidate.startswith("ENS"):
                    ids.append(candidate)
    return ids


# =====================================================================
#  Fetch CDS sequences from Ensembl
# =====================================================================


def fetch_cds_sequences(ensembl_gene: str, protein_ids: list[str]) -> dict[str, str]:
    """Fetch CDS nucleotide sequences for a set of Ensembl protein IDs.

    For each protein ID, looks up the parent transcript and then fetches
    the coding sequence (via :func:`fetch_cds_for_protein`). Querying
    ``/sequence/id/{protein_id}?type=cds`` directly returns the protein
    sequence rather than the CDS, so the parent-transcript indirection
    is required.

    Args:
        ensembl_gene: Ensembl gene ID (kept for API compatibility; currently
            unused, but retained for future homology-based fallbacks).
        protein_ids: List of Ensembl protein stable IDs from the MSA.

    Returns:
        Dict mapping protein ID to CDS nucleotide sequence string.
    """
    del ensembl_gene  # reserved for future use
    cds: dict[str, str] = {}
    for prot_id in protein_ids:
        seq = fetch_cds_for_protein(prot_id)
        if seq:
            cds[prot_id] = seq
        else:
            log.debug("No CDS for %s", prot_id)
    log.info("Fetched CDS for %d / %d proteins.", len(cds), len(protein_ids))
    return cds


def fetch_cds_for_protein(ensembl_protein: str) -> str | None:
    """Fetch the CDS nucleotide sequence for a single Ensembl protein.

    Looks up the parent transcript via ``/lookup/id`` and then fetches
    the coding sequence via ``/sequence/id``.

    Args:
        ensembl_protein: Ensembl protein stable ID (e.g. ``"ENSP00000269305"``).

    Returns:
        CDS nucleotide string, or ``None`` if unavailable.
    """
    if not ensembl_protein:
        return None

    # Get parent transcript
    url = f"{ENSEMBL_REST}/lookup/id/{ensembl_protein}"
    data = http_get_json(url)
    transcript_id = data.get("Parent") if data else None

    if not transcript_id:
        # Fallback: try direct CDS fetch on the protein ID
        url = f"{ENSEMBL_REST}/sequence/id/{ensembl_protein}?type=cds"
        data = http_get_json(url)
        return data.get("seq") if data else None

    url = f"{ENSEMBL_REST}/sequence/id/{transcript_id}?type=cds"
    data = http_get_json(url)
    return data.get("seq") if data else None


def fetch_compara_protein_msa(ensembl_gene: str, prune_taxon: int = 7742) -> dict[str, str]:
    """Fetch a Compara aligned protein MSA keyed by Ensembl protein ID.

    The Ensembl REST FASTA output of ``/genetree/member/id`` labels
    sequences with gene IDs, whereas the tree leaves use protein IDs.
    This helper walks the JSON representation of the gene tree, which
    exposes both IDs per leaf, and returns the aligned protein sequences
    keyed by the protein ID used in the tree.

    Args:
        ensembl_gene: Ensembl gene stable ID (e.g. ``"ENSG00000179344"``).
        prune_taxon: NCBI taxon ID used to restrict the gene tree
            (default 7742 = Vertebrata; 9347 = Boreoeutheria; 9443 = Primates).

    Returns:
        Ordered mapping ``{protein_id: aligned_protein_sequence}``. Empty
        on failure.
    """
    url = (
        f"{ENSEMBL_REST}/genetree/member/id/human/{ensembl_gene}"
        f"?aligned=1&sequence=protein&prune_taxon={prune_taxon}"
    )
    data = http_get_json(url)
    if not data or "tree" not in data:
        log.warning("No Compara tree returned for %s", ensembl_gene)
        return {}

    msa: dict[str, str] = {}

    def _walk(node: dict) -> None:
        if "children" in node:
            for child in node["children"]:
                _walk(child)
            return
        seq_info = node.get("sequence", {})
        mol_seq = seq_info.get("mol_seq", {})
        if isinstance(mol_seq, dict):
            aligned = mol_seq.get("seq", "")
        else:
            aligned = mol_seq or ""
        if not aligned:
            return
        # Protein accession lives under sequence.id[*]
        protein_id = ""
        for entry in seq_info.get("id", []) or []:
            acc = entry.get("accession", "") if isinstance(entry, dict) else ""
            if acc and "P" in acc:
                protein_id = acc
                break
        if not protein_id:
            return
        msa[protein_id] = aligned

    _walk(data["tree"])
    log.info("Compara protein MSA: %d sequences for %s", len(msa), ensembl_gene)
    return msa


def extract_protein_ids_from_newick(newick: str) -> list[str]:
    """Extract Ensembl protein IDs from the leaves of a Newick tree.

    Args:
        newick: Newick tree string, as returned by :func:`fetch_gene_tree_newick`.

    Returns:
        List of Ensembl protein stable IDs (``ENS...P\\d+``) found in the
        tree, preserving order.
    """
    return re.findall(r"ENS\w*P\d+", newick or "")


def fetch_translation_overlap(
    ensembl_protein: str,
    feature: str = "transcript_variation",
) -> list[dict] | None:
    """Fetch variant features overlapping a translation from Ensembl.

    Wraps the ``/overlap/translation`` REST endpoint.

    Args:
        ensembl_protein: Ensembl protein stable ID (e.g. ``"ENSP00000269305"``).
        feature: Feature type to request, e.g. ``"transcript_variation"``
            or ``"somatic_transcript_variation"``.

    Returns:
        List of variant dicts, or ``None`` on failure.
    """
    if not ensembl_protein:
        return None

    url = (
        f"{ENSEMBL_REST}/overlap/translation/{ensembl_protein}"
        f"?feature={feature};content-type=application/json"
    )
    data = http_get_json(url)
    if data and isinstance(data, list):
        return data
    return None


def fetch_pairwise_ortholog(ensembl_gene):
    log.info("Fetching pairwise ortholog alignments for %s …", ensembl_gene)
    url = (
        f"{ENSEMBL_REST}/homology/id/human/{ensembl_gene}"
        f"?type=orthologues&aligned=1&sequence=protein&format=full"
        f"&target_taxon=7742"
    )
    return http_get_json(url)


# =====================================================================
#  4. Fetch and prune gene tree
# =====================================================================


def fetch_gene_tree_newick(ensembl_gene: str, prune_taxon: int = 7742) -> str | None:
    """Fetch the Ensembl Compara gene tree in Newick format.

    Args:
        ensembl_gene: Ensembl gene stable ID.
        prune_taxon: NCBI taxon ID used to restrict the gene tree
            (default 7742 = Vertebrata; 9347 = Boreoeutheria; 9443 = Primates).

    Returns:
        Newick tree string, or ``None`` on failure.
    """
    url = (
        f"{ENSEMBL_REST}/genetree/member/id/human/{ensembl_gene}"
        f"?nh_format=simple&prune_taxon={prune_taxon}"
    )
    try:
        req = Request(
            url,
            headers={
                "Content-Type": "text/x-nh",
                "Accept": "text/x-nh",
            },
        )
        with urlopen(req, timeout=120) as resp:
            return resp.read().decode().strip()
    except Exception as e:
        log.warning("Failed to fetch gene tree: %s", e)
        return None


# =====================================================================
#  Fetch the real MSA from Ensembl Compara gene tree
# =====================================================================


def fetch_compara_msa(ensembl_gene: str, ensembl_protein: str, outpath: Path) -> str | None:
    """
    Fetch the phylogeny-aware protein MSA from the Ensembl Compara
    gene tree via the REST API genetree/member/id endpoint.

    This is the genuine multiple sequence alignment built by the
    Compara pipeline (TreeBeST), where all species in the gene tree
    are aligned simultaneously. NOT a pairwise projection.

    The endpoint supports direct FASTA output via Content-type header.

    Parameters
    ----------
    ensembl_gene : Ensembl gene ID (e.g. ENSG00000141510)
    ensembl_protein : Ensembl protein ID (fallback for lookup)
    outpath : Path to write the FASTA MSA

    Returns
    -------
    The human sequence ID in the MSA (for Jalview mapping), or None.
    """
    query_id = ensembl_gene or ensembl_protein
    if not query_id:
        log.warning("No Ensembl gene/protein ID – cannot fetch gene tree MSA.")
        return None

    log.info("Fetching Compara gene tree MSA for %s …", query_id)

    # Request the gene tree as aligned protein FASTA
    # aligned=1 → include gaps; sequence=protein → protein sequences
    # prune_taxon=7742 → Vertebrata only (keeps MSA manageable)
    url = (
        f"{ENSEMBL_REST}/genetree/member/id/human/{query_id}"
        f"?aligned=1&sequence=protein&prune_taxon=7742"
        f"&nh_format=simple"
    )

    # Request FASTA format directly
    try:
        req = Request(
            url,
            headers={
                "Content-Type": "text/x-fasta",
                "Accept": "text/x-fasta",
            },
        )
        with urlopen(req, timeout=120) as resp:
            fasta_text = resp.read().decode()
    except Exception as e:
        log.warning("Gene tree FASTA request failed: %s", e)
        log.info("Trying JSON fallback …")
        fasta_text = None

    # ── Parse FASTA response ─────────────────────────────────
    if fasta_text:
        if fasta_text.startswith(">"):
            with open(outpath, "w") as f:
                f.write(fasta_text)

            # Count sequences and find human
            human_id = None
            n_seqs = 0
            for line in fasta_text.split("\n"):
                if line.startswith(">"):
                    n_seqs += 1
                    if "homo_sapiens" in line.lower() or "human" in line.lower():
                        human_id = line[1:].split()[0]

            log.info("Compara gene tree MSA: %d sequences → %s", n_seqs, outpath)
            return human_id

    # ── JSON fallback: parse tree structure for aligned sequences ──
    log.info("Parsing gene tree JSON for aligned sequences …")
    url_json = (
        f"{ENSEMBL_REST}/genetree/member/id/human/{query_id}"
        f"?aligned=1&sequence=protein&prune_taxon=7742"
    )
    data = http_get_json(url_json)

    if not data or "tree" not in data:
        log.warning("No gene tree data returned.")
        return None

    # Recursively extract leaf nodes with aligned sequences
    records = []
    human_id = None

    def _walk_tree(node):
        """Recursively extract leaf sequences from the Compara gene tree."""
        nonlocal human_id
        if "children" in node:
            for child in node["children"]:
                _walk_tree(child)
        else:
            # Leaf node — sequence is at node["sequence"]["mol_seq"]["seq"]
            seq_data = node.get("sequence", {})
            mol_seq = seq_data.get("mol_seq", {})
            if isinstance(mol_seq, dict):
                aligned_seq = mol_seq.get("seq", "")
            elif isinstance(mol_seq, str):
                aligned_seq = mol_seq
            else:
                aligned_seq = ""
            # Fallback: try top-level seq
            if not aligned_seq:
                aligned_seq = seq_data.get("seq", "")

            species = node.get("taxonomy", {}).get("scientific_name", "")
            if not species:
                species = node.get("species", "")

            node_id = node.get("id", {})
            if isinstance(node_id, dict):
                accession = node_id.get("accession", "")
            else:
                accession = str(node_id)

            if not aligned_seq:
                return

            clean_species = species.replace(" ", "_")
            seq_id = f"{clean_species}__{accession}"[:60] if accession else clean_species[:60]

            records.append(
                SeqRecord(
                    Seq(aligned_seq),
                    id=seq_id,
                    description=species,
                )
            )

            if "homo" in species.lower() and "sapiens" in species.lower():
                human_id = seq_id

    _walk_tree(data["tree"])

    if not records:
        log.warning("No aligned sequences found in gene tree.")
        return None

    # Write FASTA
    SeqIO.write(records, outpath, "fasta")
    log.info("Compara gene tree MSA (from JSON): %d sequences → %s", len(records), outpath)

    return human_id
