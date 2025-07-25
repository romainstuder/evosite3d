#!/usr/bin/env python3

"""Move target sequence to the top of a multiple sequence alignment."""

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Tuple

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord


def find_target_sequence(
    alignment: MultipleSeqAlignment, target_id: str, exact_match: bool = False
) -> Tuple[Optional[SeqRecord], List[SeqRecord], List[SeqRecord]]:
    """Find target sequence and separate it from other records.

    Args:
        alignment: Multiple sequence alignment
        target_id: Target sequence ID or substring to match
        exact_match: If True, require exact match; if False, allow substring match

    Returns:
        Tuple of (target_record, other_records, matched_records)
    """
    target_record = None
    other_records = []
    matched_records = []

    for record in alignment:
        if exact_match:
            is_match = record.id == target_id
        else:
            is_match = target_id in record.id

        if is_match:
            if target_record is None:
                target_record = record
            else:
                matched_records.append(record)
        else:
            other_records.append(record)

    return target_record, other_records, matched_records


def move_target_sequence_to_top(
    input_file: str,
    output_file: str,
    target_id: str,
    file_format: str = "fasta",
    exact_match: bool = False,
) -> None:
    """Move target sequence to the top of a multiple sequence alignment.

    Args:
        input_file: Path to input alignment file
        output_file: Path to output alignment file
        target_id: Target sequence ID or substring to match
        file_format: Alignment file format (default: fasta)
        exact_match: If True, require exact match; if False, allow substring match

    Raises:
        FileNotFoundError: If input file doesn't exist
        ValueError: If target ID is not found or file format is invalid
    """
    try:
        # Validate input file
        input_path = Path(input_file)
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")

        if not input_path.is_file():
            raise ValueError(f"Path is not a file: {input_file}")

        # Load the alignment
        try:
            alignment = AlignIO.read(input_file, file_format)
        except Exception as e:
            raise ValueError(f"Error reading alignment file: {e}")

        if len(alignment) == 0:
            raise ValueError("Alignment file is empty")

        # Find the target sequence
        target_record, other_records, matched_records = find_target_sequence(
            alignment, target_id, exact_match
        )

        if target_record is None:
            available_ids = [record.id for record in alignment]
            raise ValueError(
                f"Target ID '{target_id}' not found in the alignment.\n"
                f"Available IDs: {', '.join(available_ids[:5])}{'...' if len(available_ids) > 5 else ''}"
            )

        # Warn about multiple matches
        if matched_records:
            print(
                f"Warning: Found {len(matched_records)} additional matches for '{target_id}':",
                file=sys.stderr,
            )
            for record in matched_records:
                print(f"  - {record.id}", file=sys.stderr)
            print(f"Using first match: {target_record.id}", file=sys.stderr)

        # Reorder the alignment
        new_alignment = MultipleSeqAlignment([target_record] + other_records + matched_records)

        # Write the output
        try:
            AlignIO.write(new_alignment, output_file, file_format)
            print(f"Saved reordered alignment to '{output_file}' with '{target_record.id}' on top.")
            print(f"Total sequences: {len(new_alignment)}")
        except Exception as e:
            raise ValueError(f"Error writing output file: {e}")

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Move a target sequence to the top of a multiple sequence alignment file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  %(prog)s input.fasta output.fasta target_seq
  %(prog)s input.fasta output.fasta target_seq --exact
  %(prog)s input.phy output.phy target_seq -f phylip
""",
    )

    parser.add_argument("input_file", help="Input MSA file (e.g., alignment.fasta)")
    parser.add_argument("output_file", help="Output MSA file")
    parser.add_argument("target_id", help="Target sequence ID or substring to match")
    parser.add_argument("-f", "--format", default="fasta", help="MSA file format (default: fasta)")
    parser.add_argument(
        "--exact",
        action="store_true",
        help="Require exact match for target ID (default: substring match)",
    )

    args = parser.parse_args()

    # Validate arguments
    if not args.target_id.strip():
        print("Error: Target ID cannot be empty", file=sys.stderr)
        sys.exit(1)

    move_target_sequence_to_top(
        args.input_file, args.output_file, args.target_id, args.format, args.exact
    )


if __name__ == "__main__":
    main()
