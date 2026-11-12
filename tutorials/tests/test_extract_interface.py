"""Tests for tutorials/.../extract_interface.py."""

import pytest

from ..protein_design.improve_protein_protein_affinity.scripts.extract_interface import (
    extract_interface,
)


def _write_interface_file(path, header_chains: list[str]) -> None:
    """Write a fake FoldX Interface_Residues file.

    The real format has 10 lines of preamble that ``read_csv(skiprows=10)``
    drops, then a tab-separated header where the columns are residue codes
    of the form ``<aa><chain><pos>`` (e.g. ``GA45``), with the *last*
    column being a trailing label that the script discards.
    """
    preamble = "\n".join(f"# preamble {i}" for i in range(10))
    columns = "\t".join(header_chains + ["junk_last"])
    body = "\t".join(["dummy"] * (len(header_chains) + 1))
    path.write_text(preamble + "\n" + columns + "\n" + body + "\n")


def test_extract_interface_returns_positions_for_chain(tmp_path) -> None:
    f = tmp_path / "interface.tsv"
    _write_interface_file(f, ["GA45", "LA47", "VB100", "RA52"])

    assert extract_interface(str(f), "A") == ["45", "47", "52"]
    assert extract_interface(str(f), "B") == ["100"]


def test_extract_interface_unknown_chain_raises(tmp_path) -> None:
    f = tmp_path / "interface.tsv"
    _write_interface_file(f, ["GA45", "LA47"])

    with pytest.raises(KeyError):
        extract_interface(str(f), "Z")


def test_extract_interface_multi_digit_positions(tmp_path) -> None:
    f = tmp_path / "interface.tsv"
    _write_interface_file(f, ["YA1234", "WA9999"])

    assert extract_interface(str(f), "A") == ["1234", "9999"]
