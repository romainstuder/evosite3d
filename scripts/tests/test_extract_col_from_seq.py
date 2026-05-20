"""Tests for scripts/extract_col_from_seq.py."""

from scripts.extract_col_from_seq import find_alignment_coordinates


def test_no_gaps() -> None:
    # Sequence positions 1..5 == alignment indices 0..4.
    aln_start, aln_stop = find_alignment_coordinates("ACGTA", 2, 4)
    # Implementation uses i-1 for start and i+1 for stop.
    assert aln_start == 0
    assert aln_stop == 4


def test_with_internal_gaps() -> None:
    # "A-C-GTA"
    # i=0 A -> j=1; i=1 '-' -> j=1; i=2 C -> j=2 (start, aln_start = 1);
    # i=3 '-' -> j stays 2 and aln_start is overwritten to 2;
    # i=4 G -> j=3; i=5 T -> j=4 (stop, aln_stop = 6, break).
    aln_start, aln_stop = find_alignment_coordinates("A-C-GTA", 2, 4)
    assert aln_start == 2
    assert aln_stop == 6


def test_start_at_position_one() -> None:
    aln_start, aln_stop = find_alignment_coordinates("ACGT", 1, 3)
    # First non-gap reaches j == 1 at i == 0, so aln_start = -1.
    assert aln_start == -1
    assert aln_stop == 3


def test_stop_beyond_sequence_returns_zero() -> None:
    aln_start, aln_stop = find_alignment_coordinates("ACGT", 1, 99)
    # Stop never matches, so default 0 is returned.
    assert aln_stop == 0
