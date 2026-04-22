# ── Taxonomy buckets (distant → close to human) ─────────────────

TAXONOMY_ORDER = [
    "Fish",
    "Amphibia",
    "Reptilia",
    "Aves",
    "Monotremata",
    "Marsupialia",
    "Laurasiatheria",
    "Rodentia",
    "Primates",
]

TAXON_MAP = {
    "Actinopterygii": "Fish",
    "Clupeocephala": "Fish",
    "Teleostei": "Fish",
    "danio_rerio": "Fish",
    "takifugu_rubripes": "Fish",
    "oryzias_latipes": "Fish",
    "gasterosteus_aculeatus": "Fish",
    "tetraodon_nigroviridis": "Fish",
    "lepisosteus_oculatus": "Fish",
    "Amphibia": "Amphibia",
    "xenopus_tropicalis": "Amphibia",
    "Sauropsida": "Reptilia",
    "Lepidosauria": "Reptilia",
    "anolis_carolinensis": "Reptilia",
    "Aves": "Aves",
    "Neognathae": "Aves",
    "gallus_gallus": "Aves",
    "taeniopygia_guttata": "Aves",
    "Monotremata": "Monotremata",
    "ornithorhynchus_anatinus": "Monotremata",
    "Marsupialia": "Marsupialia",
    "Metatheria": "Marsupialia",
    "monodelphis_domestica": "Marsupialia",
    "Laurasiatheria": "Laurasiatheria",
    "canis_lupus_familiaris": "Laurasiatheria",
    "canis_familiaris": "Laurasiatheria",
    "felis_catus": "Laurasiatheria",
    "bos_taurus": "Laurasiatheria",
    "sus_scrofa": "Laurasiatheria",
    "equus_caballus": "Laurasiatheria",
    "Cetartiodactyla": "Laurasiatheria",
    "Carnivora": "Laurasiatheria",
    "Chiroptera": "Laurasiatheria",
    "Rodentia": "Rodentia",
    "Glires": "Rodentia",
    "mus_musculus": "Rodentia",
    "rattus_norvegicus": "Rodentia",
    "cavia_porcellus": "Rodentia",
    "Primates": "Primates",
    "Simiiformes": "Primates",
    "Hominoidea": "Primates",
    "Catarrhini": "Primates",
    "pan_troglodytes": "Primates",
    "gorilla_gorilla": "Primates",
    "pongo_abelii": "Primates",
    "macaca_mulatta": "Primates",
}


def tax_classify(species: str, taxon_level: str) -> str | None:
    """Classify a species into a broad taxonomy bucket.

    Looks up the species name and taxonomy level in ``_TAXON_MAP`` to
    assign one of the groups defined in ``TAXONOMY_ORDER`` (e.g. Fish,
    Primates). Falls back to substring matching when exact keys are
    not found.

    Args:
        species: Species name (e.g. ``"mus_musculus"``).
        taxon_level: Ensembl Compara taxonomy level string.

    Returns:
        The taxonomy bucket name, or ``None`` if unclassifiable.
    """
    key = species.lower().replace(" ", "_")
    if key in TAXON_MAP:
        return TAXON_MAP[key]
    if taxon_level in TAXON_MAP:
        return TAXON_MAP[taxon_level]
    for k, v in TAXON_MAP.items():
        if k.lower() in key:
            return v
    return None


def tax_bucket_rank(b: str) -> int:
    """Return the rank of a taxonomy bucket in ``TAXONOMY_ORDER``.

    Lower ranks correspond to more distant groups (e.g. Fish = 0).
    Unknown buckets are ranked after the last defined group.

    Args:
        b: Taxonomy bucket name (e.g. ``"Primates"``).

    Returns:
        Integer rank used for sorting by evolutionary distance.
    """
    try:
        return TAXONOMY_ORDER.index(b)
    except ValueError:
        return len(TAXONOMY_ORDER)
