#!/usr/bin/env python3

"""Extract full sequences from phylogenetic analysis files."""

import argparse
import sys
from pathlib import Path
from typing import Optional


def validate_sequence(sequence: str) -> bool:
    """Validate that sequence contains only valid amino acid characters."""
    valid_chars = set("ACDEFGHIKLMNPQRSTVWY*-")
    return all(c.upper() in valid_chars for c in sequence)


def extract_sequences(
    input_file: str,
    output_file: Optional[str] = None,
    start_marker: str = "List of extant and reconstructed sequences",
    end_marker: str = "Overall accuracy of the 195 ancestral sequences",
    min_length: int = 1,
) -> None:
    """Extract full sequences from phylogenetic analysis file.

    Args:
        input_file: Path to input file
        output_file: Optional output file path (prints to stdout if None)
        start_marker: Text marker indicating start of sequence section
        end_marker: Text marker indicating end of sequence section
        min_length: Minimum sequence length to include
    """
    try:
        input_path = Path(input_file)
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")

        if not input_path.is_file():
            raise ValueError(f"Path is not a file: {input_file}")

        output_handle = open(output_file, "w") if output_file else sys.stdout

        try:
            in_sequence_section = False
            sequences_found = 0
            sequences_skipped = 0

            with open(input_file, "r") as f:
                for line_num, line in enumerate(f, 1):
                    line = line.rstrip()

                    # Check for section end marker
                    if end_marker in line:
                        in_sequence_section = False
                        continue

                    # Check for section start marker
                    if start_marker in line:
                        in_sequence_section = True
                        continue

                    # Process sequences in the target section
                    if in_sequence_section and line and not line.startswith(" "):
                        try:
                            # Clean up the line and extract sequence data
                            cleaned_line = line.replace("e #", "e_")
                            parts = cleaned_line.split()

                            if len(parts) < 2:
                                continue

                            name = parts[0]
                            sequence = "".join(parts[1:])

                            # Validate sequence
                            if len(sequence) < min_length:
                                sequences_skipped += 1
                                print(
                                    f"Warning: Skipping sequence {name} (length {len(sequence)} < {min_length})",
                                    file=sys.stderr,
                                )
                                continue

                            if not validate_sequence(sequence):
                                sequences_skipped += 1
                                print(
                                    f"Warning: Skipping sequence {name} (contains invalid characters)",
                                    file=sys.stderr,
                                )
                                continue

                            # Output sequence in FASTA format
                            print(f">{name}", file=output_handle)
                            print(sequence, file=output_handle)
                            sequences_found += 1

                        except Exception as e:
                            print(
                                f"Warning: Error processing line {line_num}: {e}", file=sys.stderr
                            )
                            continue

            # Report results
            if sequences_found == 0:
                print("Warning: No sequences found in the input file", file=sys.stderr)
            else:
                print(f"Extracted {sequences_found} sequences", file=sys.stderr)
                if sequences_skipped > 0:
                    print(f"Skipped {sequences_skipped} sequences", file=sys.stderr)

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


def main():
    parser = argparse.ArgumentParser(
        description="Extract full sequences from phylogenetic analysis files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  %(prog)s input.txt
  %(prog)s input.txt -o sequences.fasta
  %(prog)s input.txt --min-length 50
  %(prog)s input.txt --start-marker "Custom start" --end-marker "Custom end"
""",
    )

    parser.add_argument("input_file", help="Path to the input text file")
    parser.add_argument("-o", "--output", help="Output file (prints to stdout if not specified)")
    parser.add_argument(
        "--start-marker",
        default="List of extant and reconstructed sequences",
        help="Text marker indicating start of sequence section",
    )
    parser.add_argument(
        "--end-marker",
        default="Overall accuracy of the 195 ancestral sequences",
        help="Text marker indicating end of sequence section",
    )
    parser.add_argument(
        "--min-length", type=int, default=1, help="Minimum sequence length to include (default: 1)"
    )

    args = parser.parse_args()

    # Validate arguments
    if args.min_length < 1:
        print("Error: Minimum length must be at least 1", file=sys.stderr)
        sys.exit(1)

    extract_sequences(
        args.input_file, args.output, args.start_marker, args.end_marker, args.min_length
    )


if __name__ == "__main__":
    main()
