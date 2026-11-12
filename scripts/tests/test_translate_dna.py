"""Tests for scripts/translate_dna.py."""

from scripts.translate_dna import translate_fasta


def test_translate_fasta_translates_to_stop(tmp_path, capsys) -> None:
    fasta = tmp_path / "in.fasta"
    # ATG=M, GCT=A, TAA=stop -> sequence stops at TAA, only "MA" emitted.
    fasta.write_text(">gene1\nATGGCTTAA\n")

    translate_fasta(str(fasta))

    out = capsys.readouterr().out.splitlines()
    assert out == [">gene1", "MA"]


def test_translate_fasta_multiple_records(tmp_path, capsys) -> None:
    fasta = tmp_path / "in.fasta"
    fasta.write_text(">g1\nATGAAA\n>g2\nATGTGG\n")

    translate_fasta(str(fasta))

    out = capsys.readouterr().out.splitlines()
    assert out == [">g1", "MK", ">g2", "MW"]
