"""Curated reference-species panel for CodeML site-model analyses.

The panel spans ~450 My of vertebrate divergence (fishes to mammals),
deep enough for pervasive positive selection to be detectable while
small enough to keep CodeML runtime tractable on a single gene.

Both ``fetch_ensembl_msa_tree.py`` and ``run_codeml_site_models.py``
import from this module so the species list is defined in one place.
"""

from __future__ import annotations

import re

# Keys are the species prefix of an Ensembl protein stable ID — the
# letters between ``ENS`` and ``P\d+``. Human stable IDs carry no
# species letters (``ENSP\d+``) and are represented by the empty key.
# Values are ``(scientific_name, taxonomic_group)`` tuples, kept for
# logging and documentation only.
REFERENCE_SPECIES: dict[str, tuple[str, str]] = {
    # -- Mammals -------------------------------------------------------
    "": ("Homo sapiens", "Primates"),
    "PTR": ("Pan troglodytes", "Primates"),
    "GGO": ("Gorilla gorilla", "Primates"),
    "PPY": ("Pongo abelii", "Primates"),
    "MMU": ("Macaca mulatta", "Primates"),
    "CJA": ("Callithrix jacchus", "Primates"),
    "MUS": ("Mus musculus", "Rodentia"),
    "RNO": ("Rattus norvegicus", "Rodentia"),
    "OCU": ("Oryctolagus cuniculus", "Lagomorpha"),
    "CAF": ("Canis familiaris", "Carnivora"),
    "FCA": ("Felis catus", "Carnivora"),
    "BTA": ("Bos taurus", "Cetartiodactyla"),
    "SSC": ("Sus scrofa", "Cetartiodactyla"),
    "ECA": ("Equus caballus", "Perissodactyla"),
    "OAR": ("Ovis aries", "Cetartiodactyla"),
    "CHI": ("Capra hircus", "Cetartiodactyla"),
    "TUR": ("Tursiops truncatus", "Cetartiodactyla"),
    "MYL": ("Myotis lucifugus", "Chiroptera"),
    "EEU": ("Erinaceus europaeus", "Eulipotyphla"),
    "SHR": ("Sorex araneus", "Eulipotyphla"),
    "DOR": ("Dipodomys ordii", "Rodentia"),
    "CPO": ("Cavia porcellus", "Rodentia"),
    "NLE": ("Nomascus leucogenys", "Primates"),
    "TSY": ("Tarsius syrichta", "Primates"),
    "MOD": ("Monodelphis domestica", "Marsupialia"),
    "OAN": ("Ornithorhynchus anatinus", "Monotremata"),
    # -- Birds ---------------------------------------------------------
    "GAL": ("Gallus gallus", "Aves"),
    "MGA": ("Meleagris gallopavo", "Aves"),
    "TGU": ("Taeniopygia guttata", "Aves"),
    # -- Reptiles ------------------------------------------------------
    "ACA": ("Anolis carolinensis", "Reptilia"),
    "PSI": ("Pelodiscus sinensis", "Reptilia"),
    # -- Amphibians ----------------------------------------------------
    "XET": ("Xenopus tropicalis", "Amphibia"),
    # -- Fishes --------------------------------------------------------
    "LAC": ("Latimeria chalumnae", "Sarcopterygii"),
    "LOC": ("Lepisosteus oculatus", "Actinopterygii"),
    "DAR": ("Danio rerio", "Actinopterygii"),
    "AMX": ("Astyanax mexicanus", "Actinopterygii"),
    "ONI": ("Oreochromis niloticus", "Actinopterygii"),
    "GAC": ("Gasterosteus aculeatus", "Actinopterygii"),
    "ORL": ("Oryzias latipes", "Actinopterygii"),
    "TRU": ("Takifugu rubripes", "Actinopterygii"),
    "TNI": ("Tetraodon nigroviridis", "Actinopterygii"),
}


_ENSEMBL_PROTEIN_RE = re.compile(r"^ENS([A-Z]*)P\d+$")


def species_prefix(protein_id: str) -> str | None:
    """Return the Ensembl species prefix of a protein stable ID.

    Args:
        protein_id: An Ensembl protein stable ID (e.g. ``"ENSMUSP00000001234"``
            or ``"ENSP00000407332"`` for human).

    Returns:
        The letters between ``ENS`` and ``P\\d+``, empty for human, or
        ``None`` if ``protein_id`` is not a canonical Ensembl protein ID
        (e.g. the ``MGP_CAROLIEiJ_P0044190`` strain-specific form).
    """
    m = _ENSEMBL_PROTEIN_RE.match(protein_id)
    return None if m is None else m.group(1)


def is_reference_species(protein_id: str) -> bool:
    """Return ``True`` if ``protein_id`` belongs to a species in :data:`REFERENCE_SPECIES`."""
    prefix = species_prefix(protein_id)
    return prefix is not None and prefix in REFERENCE_SPECIES
