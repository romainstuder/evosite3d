"""Tests for scripts/extract_proba_from_rst.py."""

from scripts.extract_proba_from_rst import parse_positions


def test_parse_positions_no_gaps() -> None:
    # Sequence "ACGT" with shift=0; positions list maps alignment cols -> seq positions.
    # j increments after each non-gap, so j hits 0 at i==0 (before increment), then 1 after i==0, etc.
    result = parse_positions("ACGT", target_id="X", positions=[1, 3], shift=0)
    # i=0: j=0 (not in [1,3]) -> j=1
    # i=1: j=1 (in!) -> position_map[1] = 1, then j=2
    # i=2: j=2 (no) -> j=3
    # i=3: j=3 (in!) -> position_map[3] = 3, then j=4
    assert result == {1: 1, 3: 3}


def test_parse_positions_with_gaps() -> None:
    # "A-CG"
    # i=0: j=0 (no) -> j=1
    # i=1 (gap): j=1 (in!) -> position_map[1] = 1, j unchanged
    # i=2: j=1 still (in!) -> position_map[2] = 1, then j=2
    # i=3: j=2 (no) -> j=3
    result = parse_positions("A-CG", target_id="X", positions=[1], shift=0)
    assert result == {1: 1, 2: 1}


def test_parse_positions_with_shift() -> None:
    # shift just gets subtracted in the stored value.
    result = parse_positions("ACGT", target_id="X", positions=[2], shift=1)
    # i=1: j==1, not in [2], j=2
    # i=1: j==2 after increment.
    # Re-walk:
    # i=0: j=0, not in -> j=1
    # i=1: j=1, not in [2] -> j=2
    # i=2: j=2 (in!) -> position_map[2] = 2 - 1 = 1, then j=3
    # i=3: j=3, not in -> j=4
    assert result == {2: 1}


def test_parse_positions_empty_targets() -> None:
    assert parse_positions("ACGT", target_id="X", positions=[], shift=0) == {}
