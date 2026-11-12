import argparse
import logging
from itertools import combinations
from pathlib import Path
from typing import List, Tuple

import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# Constants
FOLDX_EXECUTABLE = "./foldx"
NUM_PARALLEL_PROCESSES = 7


def generate_combinations(lst: List[str], n: int) -> List[Tuple[str, ...]]:
    """Generate unique combinations of n items from a list.

    Args:
        lst: List of items to combine
        n: Number of items per combination

    Returns:
        List of unique combinations
    """
    return list(combinations(lst, n))


def filter_same_positions(variants: List[Tuple[str, ...]]) -> List[Tuple[str, ...]]:
    """Filter out variants that mutate the same position multiple times.

    Args:
        variants: List of variant tuples (e.g., ('MA10Y', 'LA20V'))

    Returns:
        Filtered list with only variants affecting different positions
    """
    filtered = []
    for variant in variants:
        # Extract position from each mutation code (e.g., 'MA10Y' -> '10')
        positions = [mutation[2:-1] for mutation in variant]
        # Keep only if all positions are unique
        if len(set(positions)) == len(positions):
            filtered.append(variant)
    return filtered


def write_mutant_commands(
    df_best: pd.DataFrame, pdb_id: str, n_mutations: int, sort_by: List[str] = None, workdir=str
) -> None:
    """Generate FoldX BuildModel commands for n-mutant variants.

    Args:
        df_best: DataFrame containing filtered mutations
        pdb_id: PDB identifier
        n_mutations: Number of mutations per variant (2 for pairs, 3 for triplets)
        sort_by: Optional list of columns to sort by before generating combinations
    """
    variant_type = "pairs" if n_mutations == 2 else "triplets"
    output_folder = Path(workdir) / f"{variant_type}_results"
    command_file = Path(workdir) / f"command_list_{variant_type}.txt"

    # Create output directory
    output_folder.mkdir(exist_ok=True)

    # Generate combinations
    if sort_by:
        mutant_codes = df_best.sort_values(sort_by)["mutant_code"].tolist()
    else:
        mutant_codes = df_best["mutant_code"].tolist()
    variant_list = generate_combinations(mutant_codes, n_mutations)
    variant_list = filter_same_positions(variant_list)

    logger.info(f"Number of variants ({variant_type}): {len(variant_list)}")

    # Write commands to file
    with command_file.open("w") as cmd_file:
        for variant in variant_list:
            # Format mutation string
            mutations_str = ",".join(variant)
            safe_name = mutations_str.replace(",", "_")

            # Create variant-specific directory
            variant_dir = output_folder / f"{pdb_id}_Repair_{safe_name}"
            variant_dir.mkdir(exist_ok=True)

            # Write mutation list file
            mutation_list_file = variant_dir / "individual_list.txt"
            mutation_list_file.write_text(f"{mutations_str};\n")

            # Build FoldX command
            command = [
                FOLDX_EXECUTABLE,
                "--command=BuildModel",
                f"--pdb={pdb_id}_Repair.pdb",
                f"--pdb-dir={workdir}",
                f"--mutant-file={mutation_list_file}",
                f"--output-dir={variant_dir}",
            ]

            cmd_file.write(" ".join(command) + "\n")

    # Log execution instructions
    logger.info(f"cat {command_file} | xargs -P {NUM_PARALLEL_PROCESSES} -I {{}} sh -c '{{}}'")


def main() -> None:
    """Generate multimutant command files for FoldX BuildModel."""
    parser = argparse.ArgumentParser(
        description="Generate multimutant commands for FoldX BuildModel",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("pdb_id", help="PDB identifier (without .pdb extension)")
    parser.add_argument(
        "dgg_threshold",
        type=float,
        help="ddG threshold for filtering mutations (mutations with ddG <= threshold will be used)",
    )
    parser.add_argument("workdir", help="Working directory")

    args = parser.parse_args()

    pdb_id = args.pdb_id
    dgg_threshold = args.dgg_threshold
    workdir = args.workdir

    # Load and filter data
    input_file = Path(workdir) / f"{pdb_id}_pssm_output.csv"
    logger.info(f"Loading PSSM data from {input_file}")

    df = pd.read_csv(input_file, sep="\t")
    df_best = df[df["ddG"] <= dgg_threshold].copy()

    logger.info(f"Found {len(df_best)} mutations with ddG <= {dgg_threshold}")
    logger.info(f"Top mutations:\n{df_best.head()}")

    # Generate pair and triplet commands
    write_mutant_commands(df_best, pdb_id, n_mutations=2, workdir=workdir)
    write_mutant_commands(
        df_best, pdb_id, n_mutations=3, sort_by=["position", "wt"], workdir=workdir
    )


if __name__ == "__main__":
    main()
