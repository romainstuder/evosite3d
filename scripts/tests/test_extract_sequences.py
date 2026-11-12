"""Tests for scripts/extract_sequences.py."""

from pathlib import Path

from scripts.extract_sequences import extract_sequences, load_gene_list


def test_load_gene_list_basic(tmp_path: Path) -> None:
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("GENE1\nGENE2\nGENE3\n")

    assert load_gene_list(str(gene_file)) == {"GENE1", "GENE2", "GENE3"}


def test_load_gene_list_skips_blanks_and_extra_columns(tmp_path: Path) -> None:
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("GENE1\textra column\n\nGENE2 ignored\n   \n")

    assert load_gene_list(str(gene_file)) == {"GENE1", "GENE2"}


def test_load_gene_list_deduplicates(tmp_path: Path) -> None:
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("GENE1\nGENE1\nGENE2\n")

    assert load_gene_list(str(gene_file)) == {"GENE1", "GENE2"}


def test_extract_sequences_returns_only_matches(tmp_path: Path) -> None:
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">GENE1\nACGT\n>GENE2\nTTTT\n>GENE3\nAAAA\n")

    records = list(extract_sequences(str(fasta), {"GENE1", "GENE3"}))

    ids = [r.id for r in records]
    assert ids == ["GENE1", "GENE3"]
    assert str(records[0].seq) == "ACGT"
    assert str(records[1].seq) == "AAAA"


def test_extract_sequences_empty_filter(tmp_path: Path) -> None:
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">GENE1\nACGT\n>GENE2\nTTTT\n")

    assert list(extract_sequences(str(fasta), set())) == []
