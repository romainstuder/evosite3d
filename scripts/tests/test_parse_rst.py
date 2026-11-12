"""Tests for scripts/parse_rst.py."""

from scripts.parse_rst import parse_rst_file


def test_parse_rst_file_emits_node_sequences(tmp_path, capsys) -> None:
    # The script splits on a 10-space delimiter, so we reproduce that exact gap.
    spacer = " " * 10
    rst = tmp_path / "rst.txt"
    rst.write_text(
        f"header line\nnode #1{spacer}A C G T\nirrelevant line\nnode #2{spacer}T T T T\n"
    )

    parse_rst_file(str(rst))

    out = capsys.readouterr().out.splitlines()
    assert out == [">node1", "ACGT", ">node2", "TTTT"]


def test_parse_rst_file_skips_short_lines(tmp_path, capsys) -> None:
    rst = tmp_path / "rst.txt"
    # No 10-space delimiter -> split returns one element -> skipped.
    rst.write_text("node only one column\n")

    parse_rst_file(str(rst))
    assert capsys.readouterr().out == ""


def test_parse_rst_file_ignores_non_node_lines(tmp_path, capsys) -> None:
    rst = tmp_path / "rst.txt"
    rst.write_text("not a node line at all\nrandom text\n")

    parse_rst_file(str(rst))
    assert capsys.readouterr().out == ""
