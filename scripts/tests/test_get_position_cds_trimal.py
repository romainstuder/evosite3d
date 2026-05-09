"""Tests for scripts/get_position_cds_trimal.py."""

import pytest

from scripts.get_position_cds_trimal import (
    build_position_mapping,
    load_trimal_columns,
    parse_site_indices,
)


def test_parse_site_indices_comma_separated() -> None:
    assert parse_site_indices("0,3,5") == [0, 3, 5]


def test_parse_site_indices_space_separated() -> None:
    assert parse_site_indices("0 3 5") == [0, 3, 5]


def test_parse_site_indices_single_value() -> None:
    assert parse_site_indices("7") == [7]


def test_parse_site_indices_invalid_raises() -> None:
    with pytest.raises(ValueError):
        parse_site_indices("0,abc")


def test_build_position_mapping_no_gaps() -> None:
    sequences = {"gene1": "ACGT"}
    mapping = build_position_mapping(sequences)
    # 1-based alignment positions -> 1-based CDS positions.
    assert mapping == {"gene1": {1: 1, 2: 2, 3: 3, 4: 4}}


def test_build_position_mapping_with_gaps() -> None:
    sequences = {"gene1": "A-CG-T"}
    mapping = build_position_mapping(sequences)
    # Gaps at align positions 2 and 5 are skipped.
    assert mapping == {"gene1": {1: 1, 3: 2, 4: 3, 6: 4}}


def test_build_position_mapping_multiple_genes() -> None:
    sequences = {"gene1": "AC", "gene2": "A-"}
    mapping = build_position_mapping(sequences)
    assert mapping == {
        "gene1": {1: 1, 2: 2},
        "gene2": {1: 1},
    }


def test_load_trimal_columns(tmp_path) -> None:
    trimal_file = tmp_path / "trimal.txt"
    trimal_file.write_text(
        "Some preamble line\n# 0, 3, 6, 9, 12\nOther content that should be ignored\n"
    )
    assert load_trimal_columns(str(trimal_file)) == [0, 3, 6, 9, 12]


def test_load_trimal_columns_uses_last_hash_line(tmp_path) -> None:
    trimal_file = tmp_path / "trimal.txt"
    trimal_file.write_text("# 1, 2, 3\n# 10, 20, 30\n")
    # Implementation overwrites on each '#' line, so the last one wins.
    assert load_trimal_columns(str(trimal_file)) == [10, 20, 30]
