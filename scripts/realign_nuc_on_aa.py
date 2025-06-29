#!/usr/bin/env python3
"""
Script to align nucleotide sequences based on a multiple sequence alignment of amino acids.

This script uses BioPython's alignment classes to properly handle MSA data.

Requirements:
    - biopython (pip install biopython)

Usage:
    python align_nucleotides.py protein_msa.fasta nucleotides.fasta output.fasta
"""

import argparse
import sys
from collections import OrderedDict
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_sequences(filename):
    """Read nucleotides sequences from FASTA file using BioPython."""
    try:
        sequences = OrderedDict()
        for record in SeqIO.parse(filename, "fasta"):
            sequences[record.id] = record
        return sequences
    except Exception as e:
        raise IOError(f"Error reading {filename}: {e}")


def read_alignment(filename):
    """Read alignment from FASTA file using BioPython AlignIO."""
    try:
        alignment = AlignIO.read(filename, "fasta")
        return alignment
    except Exception as e:
        raise IOError(f"Error reading alignment {filename}: {e}")


def write_alignment(alignment, filename):
    """Write alignment to FASTA file using BioPython AlignIO."""
    try:
        AlignIO.write(alignment, filename, "fasta")
    except Exception as e:
        raise IOError(f"Error writing alignment {filename}: {e}")


def validate_alignment_and_sequences(protein_alignment, nucleotide_records):
    """Validate that nucleotide sequences correspond to protein alignment."""
    errors = []
    warnings = []

    # Get protein sequence IDs from alignment
    protein_ids = {record.id for record in protein_alignment}
    nucleotide_ids = set(nucleotide_records.keys())

    # Check for missing nucleotide sequences
    missing_nucleotides = protein_ids - nucleotide_ids
    for prot_id in missing_nucleotides:
        errors.append(f"Missing nucleotide sequence for: {prot_id}")

    # Check for extra nucleotide sequences
    extra_nucleotides = nucleotide_ids - protein_ids
    for nt_id in extra_nucleotides:
        warnings.append(f"Nucleotide sequence {nt_id} has no corresponding protein sequence")

    # Validate sequence lengths for matching sequences
    for record in protein_alignment:
        if record.id in nucleotide_records:
            # Remove gaps from protein sequence to get original length
            protein_nogaps = str(record.seq).replace("-", "")
            nucleotide_seq = str(nucleotide_records[record.id].seq)

            # Check if nucleotide length is 3x protein length (accounting for stop codons)
            expected_nt_length = len(protein_nogaps) * 3
            actual_nt_length = len(nucleotide_seq)

            if (
                actual_nt_length != expected_nt_length
                and actual_nt_length != expected_nt_length + 3
            ):
                errors.append(
                    f"Length mismatch for {record.id}: protein={len(protein_nogaps)} AA, "
                    f"nucleotide={actual_nt_length} nt (expected {expected_nt_length} or "
                    f"{expected_nt_length + 3})"
                )

    return errors, warnings


def verify_translation(nucleotide_seq, protein_seq, genetic_code=1):
    """Verify that nucleotide sequence translates to the expected protein using BioPython."""
    try:
        # Use BioPython's translation
        nt_seq_obj = Seq(nucleotide_seq.replace("-", ""))
        translated = nt_seq_obj.translate(table=genetic_code, to_stop=False)

        # Remove stop codons and gaps for comparison
        translated_clean = str(translated).replace("*", "")
        protein_clean = protein_seq.replace("-", "").replace("*", "")

        return translated_clean == protein_clean
    except Exception as e:
        print(f"Warning: Translation verification failed: {e}")
        return False


def create_codon_alignment(protein_alignment, nucleotide_records, genetic_code=1, verbose=False):
    """Create nucleotide MultipleSeqAlignment based on protein alignment."""
    aligned_records = []

    for protein_record in protein_alignment:
        if protein_record.id not in nucleotide_records:
            if verbose:
                print(f"Warning: No nucleotide sequence found for {protein_record.id}")
            continue

        nucleotide_record = nucleotide_records[protein_record.id]
        protein_aligned = str(protein_record.seq)
        nucleotide_seq = str(nucleotide_record.seq)

        # Remove gaps from protein to get original sequence
        protein_original = protein_aligned.replace("-", "")

        # Verify translation if possible
        if len(nucleotide_seq) >= len(protein_original) * 3:
            if not verify_translation(
                nucleotide_seq[: len(protein_original) * 3], protein_original, genetic_code
            ):
                if verbose:
                    print(f"Warning: Translation mismatch for {protein_record.id}")

        # Build aligned nucleotide sequence
        aligned_nt = []
        nt_pos = 0

        for aa in protein_aligned:
            if aa == "-":
                # Insert gap (3 nucleotides worth)
                aligned_nt.append("---")
            else:
                # Add the corresponding codon
                if nt_pos + 3 <= len(nucleotide_seq):
                    codon = nucleotide_seq[nt_pos : nt_pos + 3]
                    aligned_nt.append(codon)
                    nt_pos += 3
                else:
                    # Handle case where nucleotide sequence is shorter
                    remaining = len(nucleotide_seq) - nt_pos
                    if remaining > 0:
                        partial_codon = nucleotide_seq[nt_pos:]
                        partial_codon += "N" * (3 - remaining)
                        aligned_nt.append(partial_codon)
                    else:
                        aligned_nt.append("NNN")
                    nt_pos += 3

        # Create new SeqRecord with aligned nucleotide sequence
        aligned_seq = Seq("".join(aligned_nt))
        aligned_record = SeqRecord(
            aligned_seq,
            id=nucleotide_record.id,
            description=f"Codon alignment | {nucleotide_record.description}",
        )
        aligned_records.append(aligned_record)

    # Create MultipleSeqAlignment object
    if aligned_records:
        nucleotide_alignment = MultipleSeqAlignment(aligned_records)
        return nucleotide_alignment
    else:
        return None


def get_alignment_stats(alignment):
    """Calculate alignment statistics using BioPython alignment methods."""
    if not alignment or len(alignment) == 0:
        return {}

    alignment_length = alignment.get_alignment_length()
    num_sequences = len(alignment)

    # Calculate gap statistics
    total_positions = alignment_length * num_sequences
    gap_count = sum(str(record.seq).count("-") for record in alignment)
    gap_percentage = (gap_count / total_positions) * 100 if total_positions > 0 else 0

    # Calculate conservation (positions with no gaps)
    conserved_positions = 0
    for i in range(alignment_length):
        column = alignment[:, i]  # Get column at position i
        if "-" not in column:
            conserved_positions += 1

    conservation_percentage = (
        (conserved_positions / alignment_length) * 100 if alignment_length > 0 else 0
    )

    return {
        "num_sequences": num_sequences,
        "alignment_length": alignment_length,
        "codon_length": alignment_length // 3,
        "gap_percentage": gap_percentage,
        "conservation_percentage": conservation_percentage,
    }


def print_alignment_summary(alignment, verbose=False):
    """Print a summary of the alignment using BioPython methods."""
    if not alignment:
        print("No alignment to summarize")
        return

    print("\nAlignment Summary:")
    print(f"  Sequences: {len(alignment)}")
    print(f"  Length: {alignment.get_alignment_length()} nucleotides")

    if verbose:
        print(f"  Sequence IDs: {[record.id for record in alignment]}")

        # Show first few positions of alignment
        if alignment.get_alignment_length() > 0:
            print("\nFirst 60 positions:")
            for record in alignment:
                seq_preview = str(record.seq)[:60]
                print(f"  {record.id:<15}: {seq_preview}")


def main():
    parser = argparse.ArgumentParser(
        description="Align nucleotide sequences based on protein MSA using BioPython alignment "
        "classes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script properly uses BioPython's MultipleSeqAlignment and AlignIO classes
for handling multiple sequence alignments.

Examples:
    python align_nucleotides.py proteins.fasta nucleotides.fasta aligned_nt.fasta
    python align_nucleotides.py -v -t 11 proteins.fasta nucleotides.fasta aligned_nt.fasta

Requirements:
    pip install biopython
        """,
    )

    parser.add_argument("protein_msa", help="Protein multiple sequence alignment (FASTA)")
    parser.add_argument("nucleotides", help="Nucleotide sequences (FASTA)")
    parser.add_argument("output", help="Output aligned nucleotide sequences (FASTA)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument(
        "-t",
        "--table",
        type=int,
        default=1,
        help="Genetic code table number (default: 1, standard code)",
    )
    parser.add_argument(
        "--check-files", action="store_true", help="Check if input files exist before processing"
    )
    parser.add_argument("--summary", action="store_true", help="Print detailed alignment summary")

    args = parser.parse_args()

    # Check if BioPython is available
    try:
        from Bio import __version__ as bio_version

        if args.verbose:
            print(f"Using BioPython version: {bio_version}")
    except ImportError:
        print("Error: BioPython is required. Install with: pip install biopython")
        sys.exit(1)

    # Check if files exist
    if args.check_files:
        for file_path in [args.protein_msa, args.nucleotides]:
            if not Path(file_path).exists():
                print(f"Error: File not found - {file_path}")
                sys.exit(1)

    try:
        # Read protein alignment using AlignIO
        if args.verbose:
            print(f"Reading protein MSA from: {args.protein_msa}")
        protein_alignment = read_alignment(args.protein_msa)

        # Read nucleotide sequences
        if args.verbose:
            print(f"Reading nucleotide sequences from: {args.nucleotides}")
        nucleotide_records = read_sequences(args.nucleotides)

        if args.verbose:
            print(f"Found protein alignment with {len(protein_alignment)} sequences")
            print(f"Protein alignment length: {protein_alignment.get_alignment_length()} positions")
            print(f"Found {len(nucleotide_records)} nucleotide sequences")
            print(f"Using genetic code table: {args.table}")

        # Validate sequences
        errors, warnings = validate_alignment_and_sequences(protein_alignment, nucleotide_records)

        if errors:
            print("Validation errors:")
            for error in errors:
                print(f"  {error}")
            sys.exit(1)

        if warnings and args.verbose:
            print("Warnings:")
            for warning in warnings:
                print(f"  {warning}")

        # Create nucleotide alignment
        if args.verbose:
            print("Creating nucleotide alignment based on protein MSA...")

        nucleotide_alignment = create_codon_alignment(
            protein_alignment, nucleotide_records, genetic_code=args.table, verbose=args.verbose
        )

        if not nucleotide_alignment:
            print("Error: No sequences could be aligned")
            sys.exit(1)

        # Write output using AlignIO
        if args.verbose:
            print(f"Writing aligned sequences to: {args.output}")
        write_alignment(nucleotide_alignment, args.output)

        # Print statistics
        stats = get_alignment_stats(nucleotide_alignment)
        print(f"Successfully created nucleotide alignment with {stats['num_sequences']} sequences")
        print(
            f"Alignment length: {stats['alignment_length']} nucleotides "
            f"({stats['codon_length']} codons)"
        )

        if args.verbose:
            print(f"Gap percentage: {stats['gap_percentage']:.2f}%")
            print(f"Conservation percentage: {stats['conservation_percentage']:.2f}%")

            # Validate output file
            if Path(args.output).exists():
                output_size = Path(args.output).stat().st_size
                print(f"Output file size: {output_size} bytes")

        # Print alignment summary if requested
        if args.summary:
            print_alignment_summary(nucleotide_alignment, verbose=args.verbose)

    except IOError as e:
        print(f"File I/O error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        if args.verbose:
            import traceback

            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
