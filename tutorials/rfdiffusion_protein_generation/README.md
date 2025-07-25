# RFDiffusion Tutorial: Protein Structure Generation

## Overview

RFDiffusion is a deep learning method for structure generation of proteins, protein complexes, and
other biomolecules. It uses denoising diffusion probabilistic models trained on protein structures
to generate novel protein backbones with specified constraints.

## Installation

### Prerequisites

- Python 3.10+
- CUDA-capable GPU (recommended)
- Git

### Step 1: Clone the Repository

Following instructions from there:
https://github.com/romainstuder/RFdiffusion/blob/main/README-macos.md

### Step 4: Download Model Weights

```bash
mkdir models
cd models
wget http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt
cd ..
```

## Basic Usage

### Command Structure

```bash
python scripts/run_inference.py \
    --config_path <config_file> \
    --model_dir <path_to_models> \
    --output_dir <output_directory> \
    [additional_options]
```

### Key Parameters

- `--length`: Target protein length
- `--num_designs`: Number of designs to generate
- `--pdb`: Input PDB file (for conditional generation)
- `--contigs`: Specify protein topology and constraints
- `--hotspots`: Specify important residue positions
- `--steps`: Number of diffusion steps (default: 200)

## Example 1: Unconditional Protein Generation

Generate a simple 100-residue protein:

```bash
python scripts/run_inference.py \
    inference.output_prefix=example1 \
    inference.input_pdb=null \
    'contig_map.contigs=[100-100]' \
    inference.num_designs=5 \
    inference.ckpt_override_path=./models/Base_ckpt.pt
```

**What this does:**

- Generates 5 different 100-residue protein structures
- Uses the base model checkpoint
- Outputs structures as `example1_0.pdb`, `example1_1.pdb`, etc.

## Example 2: Conditional Generation with Motif Scaffolding

Design a protein that incorporates a specific structural motif:

```bash
python scripts/run_inference.py \
    inference.output_prefix=scaffold_design \
    inference.input_pdb=input_motif.pdb \
    'contig_map.contigs=[10-40/A163-181/10-40]' \
    'ppi.hotspot_res=[A163,A165,A170,A175]' \
    inference.num_designs=10 \
    inference.ckpt_override_path=./models/Base_ckpt.pt
```

**What this does:**

- Takes an input PDB with a motif of interest
- Designs 10-40 residues before the motif (A163-181)
- Designs 10-40 residues after the motif
- Treats residues 163, 165, 170, 175 in chain A as hotspots
- Generates 10 designs

## Example 3: Protein-Protein Interface Design

Design a binder for a target protein:

```bash
python scripts/run_inference.py \
    inference.output_prefix=binder_design \
    inference.input_pdb=target_protein.pdb \
    'contig_map.contigs=[A1-150/0 70-100]' \
    'ppi.hotspot_res=[A25,A30,A45]' \
    inference.num_designs=20 \
    inference.ckpt_override_path=./models/Complex_base_ckpt.pt
```

**What this does:**

- Uses target protein (chain A, residues 1-150)
- Designs a 70-100 residue binder
- Focuses on hotspot residues for binding
- Uses the complex model for protein-protein interactions

## Working Example: Simple Helix Bundle

Let's create a simple 4-helix bundle protein:

### Step 1: Create the command

```bash
python scripts/run_inference.py \
    inference.output_prefix=helix_bundle \
    inference.input_pdb=null \
    'contig_map.contigs=[80-120]' \
    inference.num_designs=3 \
    inference.ckpt_override_path=./models/Base_ckpt.pt \
    denoiser.noise_scale_ca=1.0 \
    denoiser.noise_scale_frame=1.0
```

### Step 2: Expected Output

```
Output files:
- helix_bundle_0.pdb
- helix_bundle_1.pdb
- helix_bundle_2.pdb
- helix_bundle_0.trb (trajectory file)
- helix_bundle_1.trb
- helix_bundle_2.trb
```

### Step 3: Analyze Results

```bash
# View structure in PyMOL
pymol helix_bundle_0.pdb

# Check structure quality with basic metrics
python -c "
import biotite.structure.io as bsio
structure = bsio.load_structure('helix_bundle_0.pdb')
print(f'Number of residues: {len(structure)}')
print(f'Number of atoms: {structure.array_length()}')
"
```

## Configuration Files

RFDiffusion uses Hydra configuration files. Key config sections:

### inference.yaml

```yaml
defaults:
  - base_inference

inference:
  input_pdb: null
  output_prefix: "design"
  num_designs: 1
  ckpt_override_path: null

contig_map:
  contigs: null
  inpaint_seq: null
  length: null
```

## Tips for Success

1. **Start Simple**: Begin with unconditional generation before moving to complex constraints
2. **Model Selection**: Use `Base_ckpt.pt` for single proteins, `Complex_base_ckpt.pt` for protein
   complexes
3. **Length Constraints**: Realistic protein lengths are typically 50-500 residues
4. **Iteration**: Generate multiple designs and select the best ones
5. **Validation**: Always validate generated structures with folding prediction tools like
   AlphaFold2 or ChimeraX

## Troubleshooting

**Common Issues:**

- **CUDA out of memory**: Reduce batch size or protein length
- **Model not found**: Check model download paths
- **Invalid contigs**: Ensure contig syntax is correct: `[start-end]`
- **Poor quality structures**: Try adjusting noise scales or increasing diffusion steps

**Quality Checks:**

```bash
# Basic structure validation
python scripts/analyze_designs.py --input helix_bundle_0.pdb
```

## Next Steps

After generating structures:

1. **Structure Validation**: Use AlphaFold2, ESMFold, or other prediction tools
2. **Sequence Design**: Use ProteinMPNN or similar tools to design sequences
3. **Experimental Validation**: Express and characterize promising designs
4. **Iterative Design**: Refine based on experimental results

## Resources

- [RFDiffusion Paper](https://www.nature.com/articles/s41586-023-06415-8)
- [GitHub Repository](https://github.com/RosettaCommons/RFdiffusion)
- [Colab Tutorial](https://colab.research.google.com/drive/1bq8s3o4YxZUQTjbWVqDM4-k-j9BoZnCr)
- [IPD Software Suite](https://www.ipd.uw.edu/software/)

This tutorial provides a foundation for using RFDiffusion. Experiment with different parameters and
constraints to explore the full capabilities of this powerful protein design tool.
