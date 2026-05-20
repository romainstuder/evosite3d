"""Tests for tutorials/.../analyse_multivariants.py."""

from ..protein_design.improve_protein_protein_affinity.scripts.analyse_multivariants import (
    build_pymol_commands,
    extract_mutation_positions,
)


def test_extract_mutation_positions_pair() -> None:
    assert extract_mutation_positions("TB28Y_TB58R") == [28, 58]


def test_extract_mutation_positions_triplet() -> None:
    assert extract_mutation_positions("TB28Y_TB58R_SB104F") == [28, 58, 104]


def test_extract_mutation_positions_empty_string() -> None:
    assert extract_mutation_positions("") == []


def test_extract_mutation_positions_no_pattern_returns_empty() -> None:
    assert extract_mutation_positions("not_a_variant_name") == []


def test_extract_mutation_positions_multi_digit() -> None:
    assert extract_mutation_positions("AB1234C_DE5678F") == [1234, 5678]


def test_build_pymol_commands_contains_essentials() -> None:
    script = build_pymol_commands(
        pdb_file="model.pdb",
        mutation_positions=[28, 58, 104],
        chain_a="X",
        chain_b="Y",
        chain_a_name="receptor",
        chain_b_name="ligand",
    )

    assert "load model.pdb" in script
    # PyMOL selection format uses '+' as separator.
    assert "resi 28+58+104" in script
    # Custom chain identifiers and names get substituted.
    assert "chain X" in script
    assert "chain Y" in script
    assert "select receptor, chain X" in script
    assert "select ligand, chain Y" in script
    # Mutation summary line at the end.
    assert "positions: 28+58+104" in script


def test_build_pymol_commands_default_chains() -> None:
    script = build_pymol_commands(pdb_file="x.pdb", mutation_positions=[10])
    assert "chain A" in script
    assert "chain B" in script
    assert "select partner, chain A" in script
    assert "select target, chain B" in script
