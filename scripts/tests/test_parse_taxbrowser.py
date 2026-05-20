"""Tests for scripts/parse_taxbrowser.py."""

from scripts.parse_taxbrowser import (
    Node,
    _get_all_descendant_nodes,
    _keep_division,
    _keep_terminal,
    get_common_ancestor,
    get_genealogy,
)


def _make_taxonomy() -> dict[str, Node]:
    """Build a small synthetic taxonomy:

    1 (root)
    |
    +-- 2 (Eukaryota, superkingdom)
        |
        +-- 3 (Metazoa, kingdom)
        |   |
        |   +-- 4 (Primates, order)
        |   |   |
        |   |   +-- 5 (Homo sapiens, species, tip)
        |   |   +-- 6 (Pan troglodytes, species, tip)
        |   |
        |   +-- 7 (Trichoplax adhaerens, species, tip)
        |
        +-- 8 (Fungi, kingdom, tip)
    """
    edges = {
        "1": ("1", "no rank", "root", False),  # root, parent of itself
        "2": ("1", "superkingdom", "Eukaryota", False),
        "3": ("2", "kingdom", "Metazoa", False),
        "4": ("3", "order", "Primates", False),
        "5": ("4", "species", "Homo sapiens", True),
        "6": ("4", "species", "Pan troglodytes", True),
        "7": ("3", "species", "Trichoplax adhaerens", True),
        "8": ("2", "kingdom", "Fungi", True),
    }

    taxonomy: dict[str, Node] = {}
    for tax_id, (parent, division, name, is_tip) in edges.items():
        node = Node()
        node.tax_id = tax_id
        node.parent = parent
        node.division = division
        node.name = name
        node.is_tip = is_tip
        taxonomy[tax_id] = node

    # Wire children lists.
    for tax_id, node in taxonomy.items():
        if node.parent != tax_id and node.parent in taxonomy:
            taxonomy[node.parent].children.append(tax_id)

    return taxonomy


def test_get_genealogy_from_leaf_to_root() -> None:
    taxonomy = _make_taxonomy()
    assert get_genealogy(taxonomy, "5") == ["5", "4", "3", "2", "1"]


def test_get_genealogy_from_root_returns_self_loop() -> None:
    # Root's parent is itself; the function appends both occurrences before breaking.
    taxonomy = _make_taxonomy()
    assert get_genealogy(taxonomy, "1") == ["1", "1"]


def test_get_genealogy_unknown_taxid_returns_empty() -> None:
    taxonomy = _make_taxonomy()
    assert get_genealogy(taxonomy, "999") == []


def test_get_all_descendant_nodes() -> None:
    taxonomy = _make_taxonomy()
    descendants = _get_all_descendant_nodes(taxonomy, "3")
    assert set(descendants) == {"3", "4", "5", "6", "7"}


def test_keep_terminal_filters_to_tips() -> None:
    taxonomy = _make_taxonomy()
    descendants = _get_all_descendant_nodes(taxonomy, "3")
    tips = _keep_terminal(taxonomy, descendants)
    assert set(tips) == {"5", "6", "7"}


def test_keep_division_filters_by_division() -> None:
    taxonomy = _make_taxonomy()
    descendants = _get_all_descendant_nodes(taxonomy, "2")
    species = _keep_division(taxonomy, descendants, "species")
    assert set(species) == {"5", "6", "7"}


def test_get_common_ancestor_two_species() -> None:
    taxonomy = _make_taxonomy()
    # Homo sapiens (5) and Trichoplax (7) share Metazoa (3).
    assert get_common_ancestor(taxonomy, ["5", "7"]) == "3"


def test_get_common_ancestor_three_species() -> None:
    taxonomy = _make_taxonomy()
    # Homo sapiens (5), Pan troglodytes (6) share Primates (4),
    # add Trichoplax (7) and the answer falls back to Metazoa (3).
    assert get_common_ancestor(taxonomy, ["5", "6", "7"]) == "3"


def test_get_common_ancestor_with_self() -> None:
    taxonomy = _make_taxonomy()
    assert get_common_ancestor(taxonomy, ["5", "5"]) == "5"
