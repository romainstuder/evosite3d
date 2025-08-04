#!/usr/bin/env python3
"""
FASTA to PIR Format Converter for MODELLER

This script converts FASTA format sequences to PIR format suitable for MODELLER.
PIR format requirements:
- Header line starts with >P1; followed by sequence ID
- Second line contains sequence type and description
- Sequence lines follow, ending with '*'
"""

import argparse
import sys
from pathlib import Path


def parse_fasta(fasta_file):
    """
    Parse FASTA file and return list of (header, sequence) tuples.
    """
    sequences = []
    current_header = None
    current_sequence = []

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # Save previous sequence if exists
                if current_header is not None:
                    sequences.append((current_header, "".join(current_sequence)))
                # Start new sequence
                current_header = line[1:]  # Remove '>' character
                current_sequence = []
            elif line:
                current_sequence.append(line)

    # Don't forget the last sequence
    if current_header is not None:
        sequences.append((current_header, "".join(current_sequence)))

    return sequences


def format_pir_sequence(sequence, line_length=60):
    """
    Format sequence with specified line length and add terminating '*'.
    """
    formatted_lines = []
    for i in range(0, len(sequence), line_length):
        formatted_lines.append(sequence[i : i + line_length])

    # Add terminating '*' to the last line
    if formatted_lines:
        formatted_lines[-1] += "*"
    else:
        formatted_lines = ["*"]

    return formatted_lines


def convert_to_pir(sequences, sequence_type="sequence", description=""):
    """
    Convert FASTA sequences to PIR format for MODELLER.

    Args:
        sequences: List of (header, sequence) tuples
        sequence_type: Type of sequence ("sequence" for target, "structureX" for template)
        description: Additional description text

    Returns:
        List of PIR formatted strings
    """
    pir_entries = []

    for i, (header, sequence) in enumerate(sequences):
        # Extract sequence ID from FASTA header (first word)
        seq_id = header.split()[0]

        # Create PIR header
        pir_header = f">P1;{seq_id}"

        # Create description line
        if description:
            desc_line = f"{sequence_type}:{seq_id}:FIRST:@:LAST :@:::: {description}"
        else:
            desc_line = f"{sequence_type}:{seq_id}:FIRST:@:LAST :@::::"

        # Format sequence
        formatted_seq_lines = format_pir_sequence(sequence)

        # Combine all parts
        pir_entry = [pir_header, desc_line] + formatted_seq_lines
        pir_entries.extend(pir_entry)

        # Add blank line between entries (except for last entry)
        if i < len(sequences) - 1:
            pir_entries.append("")

    return pir_entries


def main():
    parser = argparse.ArgumentParser(
        description="Convert FASTA format to PIR format for MODELLER",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert target sequence
  python fasta_to_pir.py input.fasta -o output.pir -t sequence

  # Convert template structure
  python fasta_to_pir.py template.fasta -o template.pir -t structureX -d "Template structure"

  # Specify custom line length
  python fasta_to_pir.py input.fasta -o output.pir -l 80
        """,
    )

    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("-o", "--output", help="Output PIR file (default: stdout)")
    parser.add_argument(
        "-t",
        "--type",
        default="sequence",
        choices=["sequence", "structureX", "structureN"],
        help="Sequence type (default: sequence)",
    )
    parser.add_argument("-d", "--description", default="", help="Additional description text")
    parser.add_argument(
        "-l",
        "--line-length",
        type=int,
        default=60,
        help="Maximum line length for sequences (default: 60)",
    )

    args = parser.parse_args()

    # Check if input file exists
    if not Path(args.input).exists():
        print(f"Error: Input file '{args.input}' not found.", file=sys.stderr)
        sys.exit(1)

    try:
        # Parse FASTA file
        sequences = parse_fasta(args.input)

        if not sequences:
            print("Error: No sequences found in input file.", file=sys.stderr)
            sys.exit(1)

        print(f"Found {len(sequences)} sequence(s) in input file.", file=sys.stderr)

        # Convert to PIR format
        pir_lines = convert_to_pir(sequences, args.type, args.description)

        # Write output
        if args.output:
            with open(args.output, "w") as f:
                for line in pir_lines:
                    f.write(line + "\n")
            print(f"PIR file written to: {args.output}", file=sys.stderr)
        else:
            for line in pir_lines:
                print(line)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
