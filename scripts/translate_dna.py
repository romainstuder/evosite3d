#!/usr/bin/env python3

"""Translate DNA sequences from FASTA files to protein sequences."""

import argparse
import sys
from pathlib import Path
from typing import List, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction


def is_dna_sequence(sequence: str) -> bool:
    """Check if sequence appears to be DNA based on nucleotide content."""
    dna_chars = set("ATCGN")
    sequence_upper = sequence.upper()
    dna_content = sum(1 for char in sequence_upper if char in dna_chars)
    return dna_content / len(sequence) > 0.85 if sequence else False


def translate_sequence(
    sequence: Seq, reading_frame: int = 1, genetic_code: int = 1, to_stop: bool = True
) -> Seq:
    """Translate a DNA sequence with specified parameters.

    Args:
        sequence: DNA sequence to translate
        reading_frame: Reading frame (1, 2, 3, -1, -2, -3)
        genetic_code: Genetic code table number (default: 1 - standard code)
        to_stop: Stop translation at first stop codon

    Returns:
        Translated protein sequence
    """
    # Handle negative reading frames (reverse complement)
    if reading_frame < 0:
        sequence = sequence.reverse_complement()
        frame_offset = abs(reading_frame) - 1
    else:
        frame_offset = reading_frame - 1

    # Apply frame offset
    if frame_offset >= len(sequence):
        return Seq("")

    dna_seq = sequence[frame_offset:]

    # Translate the sequence
    try:
        protein_seq = dna_seq.translate(table=genetic_code, to_stop=to_stop)
        return protein_seq
    except Exception as e:
        raise ValueError(f"Translation error: {e}")


def translate_fasta(
    file_path: str,
    output_file: Optional[str] = None,
    reading_frames: List[int] = [1],
    genetic_code: int = 1,
    to_stop: bool = True,
    min_length: int = 10,
) -> None:
    """Translate DNA sequences from a FASTA file to protein sequences.

    Args:
        file_path: Path to input FASTA file
        output_file: Optional output file path (prints to stdout if None)
        reading_frames: List of reading frames to translate (1-6)
        genetic_code: Genetic code table number
        to_stop: Stop translation at first stop codon
        min_length: Minimum protein length to include
    """
    try:
        input_path = Path(file_path)
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {file_path}")

        if not input_path.is_file():
            raise ValueError(f"Path is not a file: {file_path}")

        output_handle = open(output_file, "w") if output_file else sys.stdout

        try:
            sequences_processed = 0
            sequences_translated = 0
            sequences_skipped = 0

            with open(file_path, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    sequences_processed += 1

                    # Check if sequence appears to be DNA
                    if not is_dna_sequence(str(record.seq)):
                        sequences_skipped += 1
                        print(
                            f"Warning: Sequence {record.id} doesn't appear to be DNA, skipping",
                            file=sys.stderr,
                        )
                        continue

                    # Translate in specified reading frames
                    for frame in reading_frames:
                        try:
                            protein_seq = translate_sequence(
                                record.seq, frame, genetic_code, to_stop
                            )

                            # Skip short proteins
                            if len(protein_seq) < min_length:
                                continue

                            # Generate header
                            frame_suffix = f"_frame{frame}" if len(reading_frames) > 1 else ""
                            header = f">{record.id}{frame_suffix}"

                            # Add sequence statistics to header
                            gc_content = gc_fraction(record.seq)
                            header += f" |gc={gc_content:.1f}%|len={len(protein_seq)}aa"

                            print(header, file=output_handle)
                            print(protein_seq, file=output_handle)
                            sequences_translated += 1

                        except Exception as e:
                            print(
                                f"Warning: Error translating {record.id} frame {frame}: {e}",
                                file=sys.stderr,
                            )
                            continue

            # Report statistics
            if sequences_processed == 0:
                print("Warning: No sequences found in input file", file=sys.stderr)
            else:
                print(f"Processed {sequences_processed} sequences", file=sys.stderr)
                print(f"Translated {sequences_translated} sequences", file=sys.stderr)
                if sequences_skipped > 0:
                    print(f"Skipped {sequences_skipped} non-DNA sequences", file=sys.stderr)

        finally:
            if output_file:
                output_handle.close()
                print(f"Results saved to: {output_file}", file=sys.stderr)

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


def parse_reading_frames(value: str) -> List[int]:
    """Parse comma-separated reading frames."""
    try:
        frames = [int(f.strip()) for f in value.split(",")]
        for frame in frames:
            if frame not in range(-3, 4) or frame == 0:
                raise ValueError(f"Invalid reading frame: {frame}")
        return frames
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"Invalid reading frames: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Translate DNA sequences from a FASTA file to protein sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  %(prog)s sequences.fasta
  %(prog)s sequences.fasta -o proteins.fasta
  %(prog)s sequences.fasta -f "1,2,3" --genetic-code 2
  %(prog)s sequences.fasta -f "1,-1" --no-stop --min-length 50

Reading frames:
  1, 2, 3: Forward frames starting at positions 0, 1, 2
  -1, -2, -3: Reverse complement frames

Genetic codes:
  1: Standard genetic code
  2: Vertebrate mitochondrial
  11: Bacterial/Archaeal
  (See NCBI genetic codes for complete list)
""",
    )

    parser.add_argument("fasta_file", help="Path to the input FASTA file containing DNA sequences")
    parser.add_argument("-o", "--output", help="Output file (prints to stdout if not specified)")
    parser.add_argument(
        "-f",
        "--frames",
        type=parse_reading_frames,
        default=[1],
        help="Comma-separated reading frames to translate (default: 1)",
    )
    parser.add_argument(
        "--genetic-code",
        type=int,
        default=1,
        help="Genetic code table number (default: 1 - standard code)",
    )
    parser.add_argument(
        "--no-stop", action="store_true", help="Don't stop translation at stop codons"
    )
    parser.add_argument(
        "--min-length", type=int, default=10, help="Minimum protein length to include (default: 10)"
    )

    args = parser.parse_args()

    # Validate arguments
    if args.genetic_code < 1 or args.genetic_code > 31:
        print("Error: Genetic code must be between 1 and 31", file=sys.stderr)
        sys.exit(1)

    if args.min_length < 1:
        print("Error: Minimum length must be at least 1", file=sys.stderr)
        sys.exit(1)

    translate_fasta(
        args.fasta_file,
        args.output,
        args.frames,
        args.genetic_code,
        not args.no_stop,
        args.min_length,
    )


if __name__ == "__main__":
    main()
