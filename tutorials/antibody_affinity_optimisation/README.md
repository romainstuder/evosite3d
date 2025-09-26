# Affinity Optimization Pipeline

This repository contains a comprehensive pipeline for protein-protein interaction affinity
optimization using FoldX molecular modeling software. The pipeline enables systematic analysis of
protein mutations to improve binding affinity through computational screening.

## Overview

The pipeline implements a complete workflow for:

1. **Structure Preparation**: Download and repair PDB structures
2. **Interface Analysis**: Identify and extract protein-protein interface residues
3. **Single Mutation Screening**: Perform comprehensive PSSM analysis
4. **Multi-mutation Design**: Generate and evaluate combinatorial mutations
5. **Variant Ranking**: Score variants based on multiple biophysical properties

## Scripts

### Structure Preparation

#### `download_pdb.py`

Downloads PDB structure files from RCSB database.

```bash
python scripts/download_pdb.py <pdb_id>
```

#### `run_foldx_repair.py`

Repairs PDB structures using FoldX RepairPDB command to optimize side chain conformations.

```bash
python scripts/run_foldx_repair.py <pdb_id>
```

### Interface Analysis

#### `run_foldx_analysecomplex.py`

Analyzes protein-protein complexes using FoldX AnalyseComplex to identify interface residues.

```bash
python scripts/run_foldx_analysecomplex.py <pdb_id> <chain1> <chain2>
```

#### `extract_interface.py`

Extracts interface residues for a specific chain and outputs PyMOL selection command.

```bash
python scripts/extract_interface.py <pdb_id> <chain_id>
```

**Output**: PyMOL selection command for visualizing interface residues

### Single Mutation Analysis

#### `generate_pssm_commands_file.py`

Generates FoldX PSSM commands for systematic single-point mutagenesis analysis.

```bash
python scripts/generate_pssm_commands_file.py <pdb_id> <other_chain> <target_chain>
```

**Output**: `command_list.txt` - Parallel execution command file

#### `parse_pssm_output.py`

Parses PSSM output files and compiles ��G values for all single mutations.

```bash
python scripts/parse_pssm_output.py <pdb_id> <target_chain>
```

**Output**: `<pdb_id>_pssm_output.csv` - Comprehensive mutation effects table

#### `plot_pssm_heatmap.py`

Generates clustered heatmap visualization of mutation effects across protein positions.

```bash
python scripts/plot_pssm_heatmap.py <pdb_id>
```

**Output**: `<pdb_id>_Repair_interface_clustermap.png` - Mutation effect heatmap

### Multi-mutation Design

#### `generate_multimutant_commands_file.py`

Generates commands for building double and triple mutant combinations from beneficial single
mutations.

```bash
python scripts/generate_multimutant_commands_file.py <pdb_id> <ddg_threshold>
```

**Output**:

- `command_list_pairs.txt` - Double mutant commands
- `command_list_triplets.txt` - Triple mutant commands

#### `analyse_multivariants.py`

Compiles results from multi-mutant FoldX calculations into a ranked output file.

```bash
python scripts/analyse_multivariants.py <pdb_id>
```

**Output**: `<pdb_id>_multivariants_output.csv` - Ranked multi-mutant results

### Variant Scoring and Selection

#### `compute_scoring.py`

Performs comprehensive scoring of variants using multiple biophysical properties including:

- Binding affinity (total energy)
- Protein stability (instability index)
- Electrostatic properties (isoelectric point)
- Hydrophobicity (GRAVY score)

```bash
python scripts/compute_scoring.py <pdb_id> <target_chain>
```

**Output**:

- `<pdb_id>_variants_ranked.csv` - Top 100 ranked variants
- `all_variant_file.fasta` - FASTA sequences of all variants

### Utilities

#### `utils.py`

Core utility functions including:

- **Amino acid mappings**: 3-letter to 1-letter code conversion
- **CDR region definitions**: CDR1, CDR2, CDR3 position ranges
- **PDB parsing**: Extract sequences from PDB structures using BioPython

## Workflow Example

Complete analysis workflow for PDB ID 9FWW:

```bash
# 1. Download and prepare structure
python scripts/download_pdb.py 9FWW
python scripts/run_foldx_repair.py 9FWW

# 2. Analyze interface
python scripts/run_foldx_analysecomplex.py 9FWW A B
python scripts/extract_interface.py 9FWW B

# 3. Single mutation screening
python scripts/generate_pssm_commands_file.py 9FWW A B
cat command_list.txt | xargs -P 7 -I {} sh -c '{}'
python scripts/parse_pssm_output.py 9FWW B
python scripts/plot_pssm_heatmap.py 9FWW

# 4. Multi-mutation design
python scripts/generate_multimutant_commands_file.py 9FWW -1.0
cat command_list_pairs.txt | xargs -P 7 -I {} sh -c '{}'
cat command_list_triplets.txt | xargs -P 7 -I {} sh -c '{}'

# 5. Analyze and rank variants
python scripts/analyse_multivariants.py 9FWW
python scripts/compute_scoring.py 9FWW B
```

## Dependencies

- **FoldX**: Molecular modeling software (requires license)
- **Python packages**:
  - pandas
  - numpy
  - matplotlib
  - seaborn
  - biopython
  - scipy

## Output Files

- `<pdb_id>_pssm_output.csv`: Single mutation ��G analysis
- `<pdb_id>_multivariants_output.csv`: Multi-mutant binding energies
- `<pdb_id>_variants_ranked.csv`: Final ranked variants with scores
- `<pdb_id>_Repair_interface_clustermap.png`: Mutation effect visualization
- `all_variant_file.fasta`: Sequences for all analyzed variants

## Scoring Methodology

Final variant scoring combines multiple factors:

- **Binding affinity** (70%): FoldX total energy
- **Stability** (10%): Protein instability index
- **Electrostatics** (10%): Isoelectric point optimization
- **Hydrophobicity** (10%): GRAVY score balance

Higher scores indicate more favorable variants for experimental validation.
