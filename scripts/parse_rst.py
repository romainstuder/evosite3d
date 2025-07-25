#!/usr/bin/env python3

"""Extract and reformat node information from an RST file."""

import argparse
import sys
from pathlib import Path
from typing import Optional, TextIO


def parse_rst_file(input_file: str, output_file: Optional[str] = None) -> None:
    """Parse RST file and extract node sequences in FASTA format.

    Args:
        input_file: Path to the RST file to parse
        output_file: Optional path to output file. If None, prints to stdout

    Raises:
        FileNotFoundError: If input file doesn't exist
        IOError: If file cannot be read or written
    """
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    output_stream: TextIO = sys.stdout
    if output_file:
        output_stream = open(output_file, "w")

    try:
        with open(input_path, "r", encoding="utf-8") as file_in:
            for line_num, line in enumerate(file_in, 1):
                line = line.rstrip()
                if line.startswith("node"):
                    parts = line.split("          ")  # Deliberate spacing used as delimiter
                    if len(parts) >= 2:
                        node_id = parts[0].replace(" ", "").replace("#", "")
                        sequence = parts[1].replace(" ", "")
                        if node_id and sequence:
                            output_stream.write(f">{node_id}\n")
                            output_stream.write(f"{sequence}\n")
                    else:
                        print(f"Warning: Malformed line {line_num}: {line}", file=sys.stderr)
    except IOError as e:
        raise IOError(f"Error reading file {input_file}: {e}") from e
    finally:
        if output_file and output_stream != sys.stdout:
            output_stream.close()


def main():
    parser = argparse.ArgumentParser(description="Extract sequences from RST file nodes.")
    parser.add_argument("input_file", help="Path to the RST file.")
    parser.add_argument("-o", "--output", help="Output file path (default: stdout)")
    args = parser.parse_args()

    try:
        parse_rst_file(args.input_file, args.output)
    except (FileNotFoundError, IOError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
