# Tutorial 1: Generating a Two-Helix Protein Structure with RFdiffusion

## Introduction

RFdiffusion is a generative model for protein structure design based on denoising diffusion
probabilistic models. It can generate diverse and realistic protein backbones from scratch or with
specific structural constraints.

## Background: Diffusion Models for Protein Design

### What is RFdiffusion?

RFdiffusion adapts the principles of image generation diffusion models (like DALL-E) to protein
structures. The model learns to reverse a diffusion process that gradually adds noise to protein
structures, enabling it to generate new structures from random noise.

### Key Concepts:

- **Denoising Process**: The model learns to remove noise from corrupted protein structures
- **SE(3) Equivariance**: Maintains proper 3D geometry and physics of proteins
- **Conditioning**: Can generate structures with specific properties (length, secondary structure,
  etc.)

## Step-by-Step Tutorial

### Step 1: Environment Setup

```bash
# Activate your conda environment
conda activate rfdiffusion

# Verify installation
python -c "import torch; print(torch.cuda.is_available())"  # Should return True
```

### Step 2: Prepare Configuration File

Create a configuration file `two_helix_config.yaml`:

```yaml
inference:
  num_designs: 10 # Generate 10 different designs
  design_startnum: 0

contigmap:
  contigs:
    - "A20-25/0 A15-20" # Two helices: 20-25 residues, loop, 15-20 residues

model_runner:
  model_name: "RFdiffusion_base"

diffusion:
  T: 50 # Number of diffusion steps

potentials:
  guide_secondary_structure: True
  secondary_structure_target: "HHHHHHHHHHHHHHHHHHHH LLLL HHHHHHHHHHHHHHHH"
  ss_bias: 10.0 # Strength of secondary structure guidance
```

### Step 3: Run RFdiffusion

```bash
# Basic command structure
python run_inference.py \
    --config two_helix_config.yaml \
    --output_dir ./output/two_helix_designs \
    --device cuda:0

# Alternative with more control
python scripts/run_inference.py \
    hydra.output_subdir=./output/two_helix_designs \
    inference.num_designs=10 \
    'contigmap.contigs=[A20-25/0 A15-20]' \
    inference.ckpt_path=./weights/RFdiffusion_base.pt \
    diffusion.T=50
```

### Step 4: Understanding the Output

RFdiffusion will generate PDB files in your output directory:

- `design_0.pdb` through `design_9.pdb`: Your generated structures
- `design_*.trb` files: Trajectory files containing metadata

### Step 5: Analyse Generated Structures

### Step 6: Visualise Results

```bash
# Open in PyMOL
pymol ./output/two_helix_designs/design_0.pdb

# PyMOL commands for visualisation
# Color by secondary structure
as cartoon
dss
color red, ss h
color yellow, ss s
color green, ss l+''

# Align all designs to see diversity
for i in range(1, 10):
    cmd.load(f'./output/two_helix_designs/design_{i}.pdb', f'design_{i}')
    cmd.align(f'design_{i}', 'design_0')
```

### Step 7: Select Best Design

Criteria for selection:

1. **Secondary Structure**: Clear two-helix architecture
2. **Compactness**: Reasonable radius of gyration
3. **Loop Geometry**: Proper connection between helices
4. **No Clashes**: Check for steric clashes

select_best_design.py

## Understanding RFdiffusion Parameters

### Key Parameters Explained:

1. **Contigs Format**: `'A20-25/0 A15-20'`

   - `A`: Chain identifier
   - `20-25`: Generate between 20-25 residues
   - `/0`: No gap (0 residues) between segments
   - Space separates different segments

2. **Diffusion Steps (T)**:

   - More steps = higher quality but slower
   - Typical range: 25-100 steps

3. **Secondary Structure Bias**:
   - `H`: Helix, `E`: Sheet, `L`: Loop
   - Higher bias = stronger enforcement

## Common Issues and Solutions

### Issue 1: Out of Memory

```bash
# Reduce batch size or number of designs
inference.num_designs=5
# Or use CPU (much slower)
--device cpu
```

### Issue 2: Poor Secondary Structure

```yaml
# Increase secondary structure bias
potentials:
  ss_bias: 20.0 # Increase from 10.0
```

### Issue 3: Unrealistic Structures

```yaml
# Add more guidance
potentials:
  guide_radius_of_gyration: True
  radius_of_gyration_target: 12.0 # Angstroms
```

## Expected Output

Your final structure should have:

- Two clear alpha helices
- Total length: 35-45 residues
- Compact fold with helices interacting
- Smooth backbone without kinks

Save your best design as `two_helix_scaffold.pdb` for use in Tutorial 2.

## Next Steps

With your scaffold structure ready, proceed to Tutorial 2 where you'll use ProteinMPNN to design an
amino acid sequence that will fold into this structure.

## Additional Resources

- RFdiffusion Paper: [Watson et al., 2023](https://www.nature.com/articles/s41586-023-06415-8)
- Official Repository: https://github.com/RosettaCommons/RFdiffusion
- Troubleshooting Guide: https://github.com/RosettaCommons/RFdiffusion/wiki
