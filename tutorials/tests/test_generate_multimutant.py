"""Tests for tutorials/.../generate_multimutant_commands_file.py."""

from ..protein_design.improve_protein_protein_affinity.scripts.generate_multimutant_commands_file import (
    filter_same_positions,
    generate_combinations,
)


def test_generate_combinations_pairs() -> None:
    result = generate_combinations(["A", "B", "C"], 2)
    assert result == [("A", "B"), ("A", "C"), ("B", "C")]


def test_generate_combinations_triples() -> None:
    result = generate_combinations(["A", "B", "C", "D"], 3)
    assert result == [
        ("A", "B", "C"),
        ("A", "B", "D"),
        ("A", "C", "D"),
        ("B", "C", "D"),
    ]


def test_generate_combinations_n_larger_than_list_returns_empty() -> None:
    assert generate_combinations(["A"], 2) == []


def test_filter_same_positions_keeps_unique_positions() -> None:
    # Variant code format: chain (1 letter) + position (digits) + new aa (1 letter).
    # Position is everything between index 1 and index -1 (so MA10Y → 'A10').
    # Note: the impl uses mutation[2:-1], i.e. drops first 2 chars and last char.
    variants = [
        ("MA10Y", "LA20V"),  # positions '10' and '20' → kept
        ("MA10Y", "LA10V"),  # both position '10' → filtered out
    ]
    assert filter_same_positions(variants) == [("MA10Y", "LA20V")]  # ty:ignore[invalid-argument-type]


def test_filter_same_positions_triplet_with_one_dup() -> None:
    variants = [
        ("MA10Y", "LA20V", "QC30P"),
        ("MA10Y", "LA20V", "QC10P"),  # last collides with first → filter
    ]
    assert filter_same_positions(variants) == [("MA10Y", "LA20V", "QC30P")]  # ty:ignore[invalid-argument-type]


def test_filter_same_positions_empty() -> None:
    assert filter_same_positions([]) == []
