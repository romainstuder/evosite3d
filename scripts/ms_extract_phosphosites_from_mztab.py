#!/usr/bin/env python3
"""
Extract phosphosite positions from mzTab files.

This script processes mzTab files to extract phosphorylation site information,
including protein accessions, peptide sequences, and absolute phosphosite positions.
"""

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

# Constants for phosphorylation UniMod IDs
# Reference: https://www.unimod.org/modifications_view.php?editid1=21
PHOSPHO_UNIMOD_IDS = {
    "00021": "Phospho",  # Generic phosphorylation
    "00035": "Oxidation",  # Sometimes co-occurs with phospho
    "00259": "Label:13C(6)15N(4)",  # Heavy labeling
}

# Default phospho IDs to search for
DEFAULT_PHOSPHO_IDS = ["00021"]

# Logging configuration
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def validate_mztab_file(file_path: str) -> None:
    """
    Validate that the input file exists and appears to be an mzTab file.

    Args:
        file_path: Path to the mzTab file

    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file doesn't appear to be a valid mzTab file
    """
    if not Path(file_path).exists():
        raise FileNotFoundError(f"mzTab file not found: {file_path}")

    # Check if file has required mzTab sections
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            content_preview = "".join([next(f, "") for _ in range(50)])

        # Look for PSM or PSH section indicators
        if not any(marker in content_preview for marker in ["PSM\t", "PSH\t"]):
            logger.warning("File may not be a standard mzTab file (no PSM/PSH sections found)")

    except Exception as e:
        raise ValueError(f"Error reading file: {e}")


def load_peptide_section(mztab_file: str) -> pd.DataFrame:
    """
    Load mzTab PSM (Peptide-Spectrum Match) section, skipping metadata lines.

    Args:
        mztab_file: Path to the mzTab file

    Returns:
        DataFrame containing PSM data

    Raises:
        ValueError: If file cannot be read or contains no PSM data
    """
    psm_lines = []

    try:
        with open(mztab_file, "r", encoding="utf-8") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line.startswith("PSM"):
                    parts = line.split("\t")
                    psm_lines.append(parts)
    except UnicodeDecodeError:
        logger.warning("UTF-8 decoding failed, trying latin-1")
        try:
            with open(mztab_file, "r", encoding="latin-1") as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if line.startswith("PSM"):
                        parts = line.split("\t")
                        psm_lines.append(parts)
        except Exception as e:
            raise ValueError(f"Error reading mzTab file with multiple encodings: {e}")
    except Exception as e:
        raise ValueError(f"Error reading mzTab file: {e}")

    if not psm_lines:
        raise ValueError("No PSM (Peptide-Spectrum Match) lines found in mzTab file")

    logger.info(f"Loaded {len(psm_lines)} PSM lines from mzTab file")
    return pd.DataFrame(psm_lines)


def get_column_indexes(header: List[str]) -> Dict[str, int]:
    """
    Get column indexes dynamically, handling mzTab format variations.

    Args:
        header: List of column names from mzTab header

    Returns:
        Dictionary mapping standard column names to their indexes

    Raises:
        ValueError: If required columns are not found
    """
    required_columns = {
        "sequence": ["sequence", "pep_sequence", "peptide_sequence"],
        "accession": ["accession", "protein_accession", "accessions"],
        "start": ["start", "pep_start", "start_position", "pep_exp_mr"],
        "modifications": ["modifications", "modification", "mods", "variable_modifications"],
    }

    column_map = {}

    for key, possible_names in required_columns.items():
        found = False
        for name in possible_names:
            if name in header:
                column_map[key] = header.index(name)
                found = True
                logger.debug(f"Found column '{key}' as '{name}' at index {header.index(name)}")
                break

        if not found:
            available_cols = [col for col in header if col.strip()][
                :10
            ]  # Show first 10 non-empty columns
            raise ValueError(
                f"Required column '{key}' not found. "
                f"Looked for: {possible_names}. "
                f"Available columns (first 10): {available_cols}"
            )

    return column_map


def extract_phosphosites(
    mod_str: str,
    peptide_start: int,
    peptide_sequence: str = "",
    unimod_ids: Optional[List[str]] = None,
) -> List[Dict]:
    """
    Extract phosphosite positions and details from modification string.

    Args:
        mod_str: Modification string from mzTab
        peptide_start: Starting position of peptide in protein sequence
        peptide_sequence: Peptide sequence (optional, for residue identification)
        unimod_ids: List of UniMod IDs to consider as phosphorylation

    Returns:
        List of dictionaries containing phosphosite information
    """
    if unimod_ids is None:
        unimod_ids = DEFAULT_PHOSPHO_IDS

    sites = []

    # Handle null/empty modifications
    if not mod_str or mod_str.lower() in ["null", "none", "", "n/a"]:
        return sites

    # Split modifications and process each one
    modifications = [mod.strip() for mod in mod_str.split(",") if mod.strip()]
    logger.debug(f"Processing modifications: {modifications}")

    for mod in modifications:
        # Try different regex patterns for modification strings
        patterns = [
            r"(\d+)\|UNIMOD:(\d+)",  # Standard format: position|UNIMOD:id
            r"(\d+)-UNIMOD:(\d+)",  # Alternative: position-UNIMOD:id
            r"(\d+)\[UNIMOD:(\d+)\]",  # Alternative: position[UNIMOD:id]
            r"(\d+)\|MOD:(\d+)",  # Legacy format: position|MOD:id
            r"(\d+)-MOD:(\d+)",  # Legacy alternative
            r"(\d+)\[MOD:(\d+)\]",  # Legacy alternative
        ]

        for pattern in patterns:
            match = re.search(pattern, mod)
            if match:
                try:
                    pos_in_pep = int(match.group(1))
                    unimod_id = match.group(2).zfill(5)  # Ensure 5-digit format

                    logger.debug(
                        f"Found modification: position={pos_in_pep}, unimod_id={unimod_id}"
                    )

                    if unimod_id in unimod_ids:
                        absolute_pos = peptide_start + pos_in_pep - 1
                        if absolute_pos > 0:  # Ensure positive position
                            # Try to determine residue if sequence is available
                            residue = ""
                            if peptide_sequence and 1 <= pos_in_pep <= len(peptide_sequence):
                                residue = peptide_sequence[pos_in_pep - 1]

                            sites.append(
                                {
                                    "absolute_position": absolute_pos,
                                    "peptide_position": pos_in_pep,
                                    "unimod_id": unimod_id,
                                    "residue": residue,
                                    "modification_name": PHOSPHO_UNIMOD_IDS.get(
                                        unimod_id, "Phospho"
                                    ),
                                }
                            )
                except ValueError as e:
                    logger.warning(f"Error parsing modification '{mod}': {e}")
                    continue
                break

    return sites


def get_header_from_mztab(mztab_path: str) -> List[str]:
    """
    Extract PSM header from mzTab file.

    Args:
        mztab_path: Path to the mzTab file

    Returns:
        List of column names from PSH (PSM header) line

    Raises:
        ValueError: If no PSH line is found
    """
    try:
        with open(mztab_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("PSH"):
                    return line.strip().split("\t")
    except UnicodeDecodeError:
        with open(mztab_path, "r", encoding="latin-1") as f:
            for line in f:
                if line.startswith("PSH"):
                    return line.strip().split("\t")

    raise ValueError("No PSH (PSM header) line found in mzTab file")


def extract_phosphosite_table(
    mztab_path: str, output_format: str = "simple", custom_unimod_ids: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Extract phosphosite information from mzTab file.

    Args:
        mztab_path: Path to the mzTab file
        output_format: 'simple' or 'detailed' output format
        custom_unimod_ids: Custom list of UniMod IDs to search for

    Returns:
        DataFrame containing extracted phosphosite information

    Raises:
        ValueError: If file processing fails
    """
    # Validate input file
    validate_mztab_file(mztab_path)

    # Use custom UniMod IDs if provided
    unimod_ids = custom_unimod_ids or DEFAULT_PHOSPHO_IDS
    logger.info(f"Searching for UniMod IDs: {unimod_ids}")

    # Load PSM section
    try:
        with open(mztab_path, "r", encoding="utf-8") as f:
            psm_lines = [line.strip() for line in f if line.startswith("PSM")]
    except UnicodeDecodeError:
        with open(mztab_path, "r", encoding="latin-1") as f:
            psm_lines = [line.strip() for line in f if line.startswith("PSM")]

    if len(psm_lines) < 1:
        raise ValueError("mzTab file must contain at least one PSM data line")

    # Get header and parse data
    header = get_header_from_mztab(mztab_path)
    data = [line.split("\t") for line in psm_lines]

    logger.info(f"Processing {len(data)} PSM entries")

    # Get column indexes
    try:
        idx = get_column_indexes(header)
        logger.info(f"Found required columns: {list(idx.keys())}")
    except ValueError as e:
        raise ValueError(f"Error parsing mzTab header: {e}")

    # Extract phosphosite data
    site_data = []
    skipped_rows = 0
    processed_rows = 0

    for row_num, row in enumerate(data, 1):
        # Ensure row has enough columns
        if len(row) <= max(idx.values()):
            skipped_rows += 1
            continue

        try:
            seq = row[idx["sequence"]] if idx["sequence"] < len(row) else ""
            acc = row[idx["accession"]] if idx["accession"] < len(row) else ""
            start = row[idx["start"]] if idx["start"] < len(row) else ""
            mods = row[idx["modifications"]] if idx["modifications"] < len(row) else ""

            # Skip rows with missing essential data
            if not seq or not acc or not start:
                skipped_rows += 1
                continue

            # Validate start position
            try:
                pep_start = int(start)
            except ValueError:
                skipped_rows += 1
                continue

            # Extract phosphosites
            sites = extract_phosphosites(mods, pep_start, seq, unimod_ids)
            processed_rows += 1

            # Add each phosphosite to results
            for site_info in sites:
                site_entry = {
                    "Protein": acc,
                    "Phosphosite_Position": site_info["absolute_position"],
                    "Peptide": seq,
                    "Modification": site_info["modification_name"],
                }

                # Add detailed information if requested
                if output_format == "detailed":
                    site_entry.update(
                        {
                            "Peptide_Start": pep_start,
                            "Position_in_Peptide": site_info["peptide_position"],
                            "Residue": site_info["residue"],
                            "UniMod_ID": site_info["unimod_id"],
                            "Modifications_String": mods,
                        }
                    )

                site_data.append(site_entry)

        except Exception as e:
            logger.warning(f"Error processing row {row_num}: {e}")
            skipped_rows += 1
            continue

    logger.info(f"Processed {processed_rows} rows, skipped {skipped_rows} rows")

    if not site_data:
        logger.warning("No phosphosites found in the input file")
        return pd.DataFrame()

    # Create DataFrame and remove duplicates
    df = pd.DataFrame(site_data)
    df = df.drop_duplicates()
    df = df.sort_values(["Protein", "Phosphosite_Position"])
    df = df.reset_index(drop=True)

    logger.info(f"Found {len(df)} unique phosphosite entries")
    return df


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Extract phosphosite positions from mzTab files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -in phospho_results.mzTab -out phosphosites.tsv
  %(prog)s -in data.mzTab -out results.csv --format detailed
  %(prog)s -in file.mzTab --unimod-ids 00021,00035
        """,
    )

    parser.add_argument(
        "-in", "--input_file", dest="input_file", required=True, help="Input mzTab file path"
    )
    parser.add_argument(
        "-out", "--output", dest="output", help="Output file path (default: phosphosite_table.tsv)"
    )
    parser.add_argument(
        "--format",
        choices=["simple", "detailed"],
        default="simple",
        help="Output format: simple or detailed (default: simple)",
    )
    parser.add_argument(
        "--unimod-ids",
        dest="unimod_ids",
        help="Comma-separated list of UniMod IDs to search for (default: 00021)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    # Configure logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Parse custom UniMod IDs if provided
    custom_unimod_ids = None
    if args.unimod_ids:
        custom_unimod_ids = [uid.strip().zfill(5) for uid in args.unimod_ids.split(",")]
        logger.info(f"Using custom UniMod IDs: {custom_unimod_ids}")

    # Set default output filename with appropriate extension
    if args.output:
        output_file = args.output
    else:
        input_path = Path(args.input_file)
        output_file = input_path.with_suffix(".phosphosites.tsv")

    try:
        logger.info(f"Processing mzTab file: {args.input_file}")
        logger.info(f"Output format: {args.format}")

        # Extract phosphosite data
        df = extract_phosphosite_table(args.input_file, args.format, custom_unimod_ids)

        if df.empty:
            logger.warning("No phosphosites found in the input file")
            return

        # Determine file format based on extension
        file_ext = Path(output_file).suffix.lower()
        if file_ext == ".tsv":
            df.to_csv(output_file, sep="\t", index=False)
        else:
            df.to_csv(output_file, index=False)

        # Print summary
        logger.info("Extraction completed successfully!")
        logger.info(f"Found {len(df)} phosphosite entries")
        logger.info(f"Unique proteins: {df['Protein'].nunique()}")
        logger.info(f"Results saved to: {output_file}")

        # Show first few results if not in verbose mode
        if not args.verbose and len(df) > 0:
            print("\nFirst 5 results:")
            print(df.head().to_string(index=False))

        # Show summary statistics
        if len(df) > 0:
            print("\nSummary Statistics:")
            print(f"Total phosphosites: {len(df)}")
            print(f"Unique proteins: {df['Protein'].nunique()}")

            if "Residue" in df.columns and not df["Residue"].empty:
                residue_counts = df["Residue"].value_counts()
                print(f"Residue distribution: {dict(residue_counts.head())}")

    except KeyboardInterrupt:
        logger.info("Process interrupted by user")
        sys.exit(0)
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback

            logger.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()
