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

# Tier definitions:
# T1: Curated reference (GRCh/T2T-level, heavily curated, gold standard)
# T2: Modern chromosome-level (long-read + Hi-C; VGP-class)
# T3: Legacy chromosome-level (older assemblies, decent but dated)
# T4: Scaffold-level (moderate contiguity)
# T5: Fragmented / low-quality


TIER_1_REFERENCE = {
    # -- Primates (T2T / GRC deep-curated) ----------------------------
    "": ("Homo sapiens", "Primates"),
    # -- Primates (chromosome-level, deeply annotated) ----------------
    "MMU": ("Macaca mulatta", "Primates"),
    # -- Rodentia (GRC deep-curated) ----------------------------------
    "MUS": ("Mus musculus", "Rodentia"),
    "RNO": ("Rattus norvegicus", "Rodentia"),
    # -- Birds --------------------------------------------------------
    "GAL": ("Gallus gallus", "Aves"),
    # -- Fishes -------------------------------------------------------
    "DAR": ("Danio rerio", "Actinopterygii"),
}

TIER_2_CHROMOSOME = {
    # -- Primates -----------------------------------------------------
    "PTR": ("Pan troglodytes", "Primates"),
    "GGO": ("Gorilla gorilla", "Primates"),
    "PPY": ("Pongo abelii", "Primates"),
    "PPA": ("Pan paniscus", "Primates"),
    "NLE": ("Nomascus leucogenys", "Primates"),
    "CJA": ("Callithrix jacchus", "Primates"),
    # -- Rodentia -----------------------------------------------------
    "CPO": ("Cavia porcellus", "Rodentia"),
    "ITC": ("Ictidomys tridecemlineatus", "Rodentia"),
    # -- Lagomorpha ---------------------------------------------------
    "OCU": ("Oryctolagus cuniculus", "Lagomorpha"),
    # -- Cetartiodactyla ----------------------------------------------
    "BTA": ("Bos taurus", "Cetartiodactyla"),
    "SSC": ("Sus scrofa", "Cetartiodactyla"),
    "OAR": ("Ovis aries", "Cetartiodactyla"),
    "CHI": ("Capra hircus", "Cetartiodactyla"),
    "VPA": ("Vicugna pacos", "Cetartiodactyla"),
    # -- Perissodactyla -----------------------------------------------
    "ECA": ("Equus caballus", "Perissodactyla"),
    # -- Carnivora ----------------------------------------------------
    "CAF": ("Canis familiaris", "Carnivora"),  # ROS_Cfam_1.0: scaffold N50 64 Mb
    "FCA": ("Felis catus", "Carnivora"),
    "MPU": ("Mustela putorius furo", "Carnivora"),
    "AML": ("Ailuropoda melanoleuca", "Carnivora"),
    # -- Chiroptera ---------------------------------------------------
    "PTL": ("Pteropus vampyrus", "Chiroptera"),
    "MDA": ("Myotis davidii", "Chiroptera"),
    # -- Marsupialia --------------------------------------------------
    "MOD": ("Monodelphis domestica", "Marsupialia"),
    # -- Birds --------------------------------------------------------
    "MGA": ("Meleagris gallopavo", "Aves"),
    "TGU": ("Taeniopygia guttata", "Aves"),
    # -- Fishes -------------------------------------------------------
    "ORL": ("Oryzias latipes", "Actinopterygii"),
    "GAC": ("Gasterosteus aculeatus", "Actinopterygii"),
    "ONI": ("Oreochromis niloticus", "Actinopterygii"),
    "TRU": ("Takifugu rubripes", "Actinopterygii"),
    "TNI": ("Tetraodon nigroviridis", "Actinopterygii"),
}

TIER_3_MODERATE = {
    # -- Primates -----------------------------------------------------
    "MML": ("Microcebus murinus", "Primates"),
    "OGA": ("Otolemur garnettii", "Primates"),
    # -- Scandentia ---------------------------------------------------
    "TBE": ("Tupaia belangeri", "Scandentia"),
    # -- Chiroptera ---------------------------------------------------
    "MYL": ("Myotis lucifugus", "Chiroptera"),  # scaffold N50 4.3 Mb — moved UP from Tier 4
    # -- Afrotheria ---------------------------------------------------
    "LAF": ("Loxodonta africana", "Afrotheria"),
    # -- Monotremata --------------------------------------------------
    "OAN": ("Ornithorhynchus anatinus", "Monotremata"),
    # -- Reptiles -----------------------------------------------------
    "ACA": ("Anolis carolinensis", "Reptilia"),
    "PSI": ("Pelodiscus sinensis", "Reptilia"),
    # -- Amphibians ---------------------------------------------------
    "XET": ("Xenopus tropicalis", "Amphibia"),
    # -- Fishes -------------------------------------------------------
    "LOC": ("Lepisosteus oculatus", "Actinopterygii"),
    "AMX": ("Astyanax mexicanus", "Actinopterygii"),
    "LAC": ("Latimeria chalumnae", "Sarcopterygii"),
}

TIER_4_FRAGMENTED = {
    # -- Primates -----------------------------------------------------
    "TSY": ("Tarsius syrichta", "Primates"),  # scaffold N50 ~58 kb
    # -- Rodentia -----------------------------------------------------
    "DOR": ("Dipodomys ordii", "Rodentia"),
    # -- Cetartiodactyla ----------------------------------------------
    "TUR": (
        "Tursiops truncatus",
        "Cetartiodactyla",
    ),  # scaffold N50 157 kb, 2.59× — gene-scaffolds required
    # -- Eulipotyphla -------------------------------------------------
    "EEU": ("Erinaceus europaeus", "Eulipotyphla"),  # 1.86× coverage, no N50 stated
    "SHR": ("Sorex araneus", "Eulipotyphla"),  # scaffold N50 18.7 kb — worst in list
}

REFERENCE_SPECIES = TIER_1_REFERENCE | TIER_2_CHROMOSOME


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
