"""Tests for scripts/remove_empty_col.py."""

from scripts.remove_empty_col import remove_empty_columns


def test_drops_columns_that_are_all_dashes() -> None:
    sequences = {
        "A": "AC-GT",
        "B": "AC-GT",
        "C": "AC-GT",
    }
    assert remove_empty_columns(sequences) == {
        "A": list("ACGT"),
        "B": list("ACGT"),
        "C": list("ACGT"),
    }


def test_drops_columns_that_are_all_X() -> None:
    sequences = {
        "A": "AXG",
        "B": "AXG",
    }
    assert remove_empty_columns(sequences) == {
        "A": list("AG"),
        "B": list("AG"),
    }


def test_keeps_column_with_mixed_chars() -> None:
    sequences = {
        "A": "A-G",
        "B": "AAG",
    }
    # Middle column has both '-' and 'A', so it must be kept.
    assert remove_empty_columns(sequences) == {
        "A": list("A-G"),
        "B": list("AAG"),
    }


def test_keeps_column_with_only_residues() -> None:
    sequences = {
        "A": "ACG",
        "B": "ACG",
    }
    assert remove_empty_columns(sequences) == {
        "A": list("ACG"),
        "B": list("ACG"),
    }


def test_drops_only_dash_or_X_columns_mixed() -> None:
    sequences = {
        "A": "A-XG",
        "B": "A-XG",
    }
    assert remove_empty_columns(sequences) == {
        "A": list("AG"),
        "B": list("AG"),
    }
