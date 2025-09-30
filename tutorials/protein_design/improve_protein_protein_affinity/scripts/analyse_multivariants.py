import argparse
import logging
import re
from pathlib import Path
from typing import List

import pandas as pd

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def extract_mutation_positions(variant_name: str) -> List[int]:
    """Extract mutation positions from variant name.

    Args:
        variant_name: Variant name like "TB28Y_TB58R_SB104F"

    Returns:
        List of mutation positions (e.g., [28, 58, 104])
    """
    # Pattern: one letter (chain) + digits (position) + one letter (new amino acid)
    pattern = r"[A-Z]([0-9]+)[A-Z]"
    positions = re.findall(pattern, variant_name)
    return [int(pos) for pos in positions]


def build_pymol_commands(
    pdb_file: str,
    mutation_positions: List[int],
    chain_a: str = "A",
    chain_b: str = "B",
    chain_a_name: str = "partner",
    chain_b_name: str = "target",
) -> str:
    """
    Generate PyMOL visualisation commands for a variant.

    Args:
        pdb_file: Path to PDB file
        mutation_positions: List of mutated residue positions
        chain_a: Chain A identifier
        chain_b: Chain B identifier (chain with mutations)
        chain_a_name: Display name for chain A
        chain_b_name: Display name for chain B

    Returns:
        PyMOL script as string
    """
    # Format mutation positions for PyMOL selection
    positions_str = "+".join(str(pos) for pos in mutation_positions)

    # TODO: Use generic elements instead of CDR.
    # TODO: Use better colours.
    command_block = f"""# PyMOL visualisation script for variant
# Generated automatically by analyse_multivariants.py

# Load structure
load {pdb_file}

# Define chains
select {chain_a_name}, chain {chain_a}
select {chain_b_name}, chain {chain_b}

# Basic visualisation
hide everything
show cartoon, {chain_a_name}
show cartoon, {chain_b_name}
color grey60, {chain_a_name}
color white, {chain_b_name}

# Define CDR regions (if applicable)
select CDR1, chain {chain_b} and resi 26-35
select CDR2, chain {chain_b} and resi 50-66
select CDR3, chain {chain_b} and resi 97-112
color cyan, CDR1
color magenta, CDR2
color olive, CDR3

# Highlight mutated positions
select mutated, chain {chain_b} and resi {positions_str}
show spheres, mutated
# color yellow, mutated

# Show disulfide bridges
select disulfides, resn CYS
show sticks, disulfides
color orange, disulfides

# Show contact residues on partner chain
select contact_residues, chain {chain_a} within 5.0 of mutated
show spheres, contact_residues
# color salmon, contact_residues

# Apply element colors to spheres
util.cnc mutated
util.cnc contact_residues

# Clean up selections
select none

# Center view on mutated residues
zoom mutated, 8

# Print summary
print "Variant loaded with mutations at positions: {positions_str}"
"""
    return command_block


def main():
    parser = argparse.ArgumentParser(
        description="Analyse multivariants from FoldX output and "
        "generate PyMOL visualisation scripts"
    )
    parser.add_argument("pdb_id", help="PDB identifier (without .pdb extension)")
    parser.add_argument(
        "--chain-a", default="A", help="Chain identifier for binding partner (default: A)"
    )
    parser.add_argument(
        "--chain-b",
        default="B",
        help="Chain identifier for target chain with mutations (default: B)",
    )
    parser.add_argument(
        "--chain-a-name",
        default="partner",
        help="Display name for chain A in PyMOL (default: partner)",
    )
    parser.add_argument(
        "--chain-b-name",
        default="target",
        help="Display name for chain B in PyMOL (default: target)",
    )
    parser.add_argument(
        "--no-pymol", action="store_true", help="Skip PyMOL visualisation script generation"
    )
    parser.add_argument("workdir", help="Working directory")

    args = parser.parse_args()

    pdb_id = args.pdb_id
    workdir = args.workdir

    df_list = []
    pymol_scripts_created = 0

    for variants_type in ["pairs", "triplets"]:
        # Load variants_type
        file = open(f"{workdir}/command_list_{variants_type}.txt", "r")
        line = file.readline()
        folder_list = []
        while line:
            line = file.readline().rstrip()
            if line != "":
                output_folder = line.split("--output-dir=")[-1]
                folder_list.append(output_folder)
        file.close()
        logger.info(f"Number of {variants_type} variants: {len(folder_list)}")

        for output_folder in folder_list:
            average_file = f"{output_folder}/Average_{pdb_id}_Repair.fxout"
            df = pd.read_csv(average_file, sep="\t", skiprows=8)

            # Extract variant name
            variant_name = output_folder.replace(
                f"{workdir}/{variants_type}_results/{pdb_id}_Repair_", ""
            )
            df["variant"] = variant_name
            df_list.append(df)

            # Generate PyMOL visualisation script (unless --no-pymol is specified)
            if not args.no_pymol:
                try:
                    # Extract mutation positions from variant name
                    mutation_positions = extract_mutation_positions(variant_name)

                    if mutation_positions:
                        # Find PDB file in the output folder
                        folder_path = Path(output_folder)
                        pdb_files = list(folder_path.glob(f"{pdb_id}_Repair_*.pdb"))

                        if pdb_files:
                            # Use the first PDB file found (usually _1.pdb)
                            pdb_file = pdb_files[0].name

                            # Generate PyMOL commands
                            pymol_script = build_pymol_commands(
                                pdb_file=pdb_file,
                                mutation_positions=mutation_positions,
                                chain_a=args.chain_a,
                                chain_b=args.chain_b,
                                chain_a_name=args.chain_a_name,
                                chain_b_name=args.chain_b_name,
                            )

                            # Write PyMOL script to variant folder
                            pymol_file = folder_path / "visualise.pml"
                            pymol_file.write_text(pymol_script)
                            pymol_scripts_created += 1

                            logger.debug(f"Created PyMOL script for {variant_name} at {pymol_file}")
                        else:
                            logger.warning(f"No PDB file found in {output_folder}")
                    else:
                        logger.warning(
                            f"Could not extract positions from variant name: {variant_name}"
                        )

                except Exception as e:
                    logger.error(f"Error generating PyMOL script for {variant_name}: {e}")

    if not args.no_pymol and pymol_scripts_created > 0:
        logger.info(f"Created {pymol_scripts_created} PyMOL visualisation scripts")

    # Combine and rank all variants
    df_output = pd.concat(df_list)
    df_output = df_output.sort_values(["total energy"])
    cols = ["variant"] + [col for col in df_output.columns if col != "variant"]
    df_output = df_output.reindex(columns=cols)

    output_file = f"{workdir}/{pdb_id}_multivariants_output.csv"
    df_output.to_csv(output_file, index=False, sep="\t")
    logger.info(f"Results written to {output_file}")

    # Display top variants
    logger.info("\n=== Top 10 Variants by Total Energy ===")
    top_variants = df_output.head(10)
    for idx, (_, row) in enumerate(top_variants.iterrows(), 1):
        variant_name = row["variant"]
        total_energy = row["total energy"]
        logger.info(f"{idx:2d}. {variant_name:40s} Total Energy: {total_energy:8.2f}")

    # Show PyMOL usage instructions
    if not args.no_pymol and pymol_scripts_created > 0:
        logger.info("\n=== PyMOL Visualization ===")
        logger.info("To visualize a variant, navigate to its folder and run:")
        logger.info("  cd <variant_folder>")
        logger.info("  pymol visualise.pml")
        logger.info("\nExample:")
        if len(df_output) > 0:
            top_variant = df_output.iloc[0]["variant"]
            # Determine folder type based on number of underscores in variant name
            folder_type = "pairs_results" if top_variant.count("_") == 1 else "triplets_results"
            logger.info(f"  cd {workdir}/{folder_type}/{pdb_id}_Repair_{top_variant}")
            logger.info("  pymol visualise.pml")


if __name__ == "__main__":
    main()
