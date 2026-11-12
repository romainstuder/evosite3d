"""Tests for scripts/extract_full_sequences.py."""

from scripts.extract_full_sequences import extract_sequences


def test_extracts_sequences_within_section(tmp_path, capsys) -> None:
    text = (
        "preamble\n"
        "List of extant and reconstructed sequences\n"
        "node #1   A C G T\n"
        "node #2   T T T T\n"
        "Overall accuracy of the 195 ancestral sequences: foo\n"
        "node #3   X X X X\n"  # after the end marker, must be skipped
    )
    inp = tmp_path / "rst.txt"
    inp.write_text(text)

    extract_sequences(str(inp))

    out = capsys.readouterr().out.splitlines()
    assert out == [">node_1", "ACGT", ">node_2", "TTTT"]


def test_skips_indented_lines(tmp_path, capsys) -> None:
    text = (
        "List of extant and reconstructed sequences\n"
        "  indented line should be skipped\n"
        "node #1   A C G T\n"
    )
    inp = tmp_path / "rst.txt"
    inp.write_text(text)

    extract_sequences(str(inp))

    out = capsys.readouterr().out.splitlines()
    assert out == [">node_1", "ACGT"]


def test_no_section_markers_yields_nothing(tmp_path, capsys) -> None:
    inp = tmp_path / "rst.txt"
    inp.write_text("just\nrandom\nlines\n")

    extract_sequences(str(inp))
    assert capsys.readouterr().out == ""
