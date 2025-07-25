#!/usr/bin/env python3

"""Extract transition between branches from phylogenetic analysis files"""

import argparse
import sys
from pathlib import Path
from typing import List, Optional


def extract_transitions(
    file_path: str,
    target_branches: List[str],
    target_positions: List[str],
    output_file: Optional[str] = None,
) -> None:
    """Extract transitions between branches from phylogenetic analysis file.

    Args:
        file_path: Path to input file
        target_branches: List of branch identifiers to extract
        target_positions: List of position identifiers to extract
        output_file: Optional output file path (prints to stdout if None)
    """
    try:
        input_path = Path(file_path)
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {file_path}")

        if not input_path.is_file():
            raise ValueError(f"Path is not a file: {file_path}")

        output_handle = open(output_file, "w") if output_file else sys.stdout

        try:
            parsing_section = False
            current_branch = ""
            matches_found = 0

            with open(file_path, "r") as file_in:
                for line_num, line in enumerate(file_in, 1):
                    line = line.rstrip()

                    # Check for section boundaries
                    if "List of extant and reconstructed sequences" in line:
                        parsing_section = False
                        continue

                    if "Summary of changes along branches" in line:
                        parsing_section = True
                        continue

                    # Process lines in the target section
                    if parsing_section:
                        tokens = line.split()
                        if not tokens:
                            continue

                        # Check if this line defines a branch
                        if tokens[0] == "Branch" and len(tokens) > 2:
                            current_branch = tokens[2]
                            continue

                        # Check if this line contains a position of interest
                        position = tokens[0]
                        if position in target_positions and current_branch in target_branches:
                            print(f"{current_branch}\t{position}\t{line}", file=output_handle)
                            matches_found += 1

            if matches_found == 0:
                print(
                    "Warning: No matches found for the specified branches and positions",
                    file=sys.stderr,
                )
            else:
                print(f"Found {matches_found} matches", file=sys.stderr)

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


def parse_list_argument(value: str) -> List[str]:
    """Parse comma-separated list argument."""
    return [item.strip() for item in value.split(",") if item.strip()]


def main():
    # Default values
    default_branches = ["208..209", "209..210", "209..265", "208..230"]
    default_positions = ["198", "204", "222", "234", "251", "270", "272", "307", "318"]

    parser = argparse.ArgumentParser(
        description="Extract transitions between branches from phylogenetic analysis files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  %(prog)s input.txt
  %(prog)s input.txt -b "208..209,209..210" -p "198,204,222"
  %(prog)s input.txt -o output.txt
""",
    )

    parser.add_argument("input_file", help="Input file to parse")
    parser.add_argument(
        "-b",
        "--branches",
        type=parse_list_argument,
        default=default_branches,
        help=f"Comma-separated list of target branches (default: {','.join(default_branches)})",
    )
    parser.add_argument(
        "-p",
        "--positions",
        type=parse_list_argument,
        default=default_positions,
        help=f"Comma-separated list of target positions (default: {','.join(default_positions)})",
    )
    parser.add_argument("-o", "--output", help="Output file (prints to stdout if not specified)")

    args = parser.parse_args()

    # Validate arguments
    if not args.branches:
        print("Error: At least one branch must be specified", file=sys.stderr)
        sys.exit(1)

    if not args.positions:
        print("Error: At least one position must be specified", file=sys.stderr)
        sys.exit(1)

    extract_transitions(args.input_file, args.branches, args.positions, args.output)


if __name__ == "__main__":
    main()
