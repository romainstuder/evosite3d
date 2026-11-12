"""Tests for tutorials/.../improve_protein_protein_affinity/scripts/utils.py."""

from ..protein_design.improve_protein_protein_affinity.scripts.utils import (
    AA_LIST,
    CDR1,
    CDR2,
    CDR3,
    D3TO1,
    AminoAcidConstants,
    AntibodyConstants,
)


def test_d3to1_has_all_20_amino_acids() -> None:
    assert len(D3TO1) == 20
    # Every value is a single uppercase letter.
    for one_letter in D3TO1.values():
        assert len(one_letter) == 1
        assert one_letter.isupper()


def test_d3to1_specific_mappings() -> None:
    assert D3TO1["ALA"] == "A"
    assert D3TO1["TRP"] == "W"
    assert D3TO1["GLY"] == "G"


def test_aa_list_length_matches_d3to1() -> None:
    assert len(AA_LIST) == len(D3TO1)


def test_aa_list_uses_three_letter_codes() -> None:
    for code in AA_LIST:
        assert len(code) == 3


def test_legacy_aliases_mirror_class_constants() -> None:
    assert D3TO1 is AminoAcidConstants.D3TO1
    assert AA_LIST is AminoAcidConstants.AA_LIST
    assert CDR1 is AntibodyConstants.CDR1
    assert CDR2 is AntibodyConstants.CDR2
    assert CDR3 is AntibodyConstants.CDR3


def test_cdr_ranges_are_kabat_ranges() -> None:
    # Kabat numbering ranges for CDR loops.
    assert list(CDR1) == list(range(26, 34))
    assert list(CDR2) == list(range(51, 59))
    assert list(CDR3) == list(range(97, 113))


def test_cdrs_do_not_overlap() -> None:
    assert set(CDR1).isdisjoint(set(CDR2))
    assert set(CDR2).isdisjoint(set(CDR3))
    assert set(CDR1).isdisjoint(set(CDR3))
