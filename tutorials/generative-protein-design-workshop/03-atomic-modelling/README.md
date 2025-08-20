# Tutorial 3: Building Full Atomic Models with MODELLER

## Introduction

MODELLER is a program for comparative protein structure modeling. In this tutorial, we'll use it in
a creative way: building a full atomic model by treating our RFdiffusion scaffold as a "template"
and threading our ProteinMPNN-designed sequence onto it.

## Background: Homology Modeling Principles

### What is MODELLER?

MODELLER constructs 3D models of proteins by:

- **Spatial Restraints**: Derives restraints from template structures
- **Optimization**: Satisfies spatial restraints while maintaining proper stereochemistry
- **Side Chain Modeling**: Places side chains based on rotamer libraries

### Our Approach:

We're using MODELLER unconventionally:

- **Template**: Our RFdiffusion backbone (Cα trace)
- **Target**: Our ProteinMPNN sequence
- **Goal**: Build full atomic coordinates

## Step-by-Step Tutorial

### Environment Setup

Activate your conda environment and install Modeller with

```shell
conda activate rfdiffusion
conda install -c salilab modeller
```

```shell
export EVOSITE3D_SCRIPTS=./utils
```

### Step 1: Prepare Input Files

First, we need to prepare our scaffold and sequence:

```shell
python $EVOSITE3D_SCRIPTS/prepare_inputs.py \
--scaffold-pdb scaffold_clean.pdb \
--sequence-file designed_sequence.fasta \
--verbose
```

### Step 2: Create MODELLER Script

build_model.py

### Step 3: Run MODELLER

```bash
# Make sure MODELLER is properly licensed

# Run the complete script
python $EVOSITE3D_SCRIPTS/build_model.py

# The script will output something like:
# Prepared files for MODELLER:
# - Template: template.pdb (45 residues)
# - Target sequence: 45 residues
# - Alignment: alignment.ali
#
# Building models...
# Model: target.B99990001.pdb
#   DOPE score: -4782.23
#   GA341 score: 0.9823
# ...
# Best model: target.B99990001.pdb (DOPE: -4782.23)
# Best model saved as: best_model.pdb
```

### Step 4: Analyse Generated Models

After MODELLER completes, you will have several model files. Let's analyze them:

```shell
python $EVOSITE3D_SCRIPTS/analyse_models.py
```

Choose the best model and save as `best_model.pdb`

```shell
cp target.B99990005.pdb best_model.pdb
```

### Step 5: Validate the Best Model

Now let's thoroughly validate our best model:

```shell
python $EVOSITE3D_SCRIPTS/validate_best_model.py
```

### Step 6: Visualise the Final Model

Use PyMOL to inspect your final model:

```bash
# Open PyMOL with both the scaffold and final model
pymol best_model.pdb scaffold_clean.pdb

# In PyMOL, execute these commands:
# 1. Align the structures
align best_model and name CA, scaffold_clean and name CA

# 2. Color by structure
hide all
show cartoon
color paleyellow, best_model and ss h
color yellow, best_model and ss s
color grey60, best_model and ss l+''

# 3. Show side chains
show sticks, best_model
hide sticks, name c+n+o

# 4. Display hydrophobic core
select hydrophobic, best_model and resn ILE+LEU+VAL+PHE+TRP+TYR+MET
show surface, hydrophobic
color orange, hydrophobic
set transparency, 0.9

# 5. Save the session
save final_design_session.pse
```

### Step 7: Export Final Results

Create a summary of your design:

```shell
python $EVOSITE3D_SCRIPTS/export_results.py
```

## Final Checklist

Before considering your design complete, verify:

- [ ] **Files exist:**

  - [ ] `two_helix_scaffold.pdb` (from Tutorial 1)
  - [ ] `designed_sequence.fasta` (from Tutorial 2)
  - [ ] `best_model.pdb` (from Tutorial 3)

- [ ] **Model quality:**

  - [ ] All atoms present (check with validation script)
  - [ ] No severe clashes (< 5)
  - [ ] Correct number of residues
  - [ ] Helical structure preserved

- [ ] **Documentation:**
  - [ ] Save PyMOL session
  - [ ] Export summary report
  - [ ] Keep all intermediate files

## Troubleshooting Common Issues

### "No models generated"

- Check that sequence length matches scaffold length
- Verify alignment file format
- Ensure MODELLER license is valid

### "Poor model quality"

- Generate more models (increase `a.ending_model`)
- Try different refinement levels
- Check input file quality

### "Missing atoms in model"

- This is usually fine for terminal residues
- For other residues, check the alignment

## Understanding MODELLER Parameters

### Key Concepts:

1. **Spatial Restraints**: Constraints that guide modeling

   - Distance restraints
   - Dihedral angle restraints
   - Secondary structure restraints

2. **Optimization Levels**:

   - `refine.very_fast`: Quick but rough
   - `refine.fast`: Balanced
   - `refine.slow`: Thorough
   - `refine.slow_large`: Most thorough

3. **Assessment Methods**:
   - **DOPE**: Statistical potential (lower is better)
   - **GA341**: Model quality (0-1, higher is better)
   - **Normalized DOPE**: Z-score (< -1 is good)

## Expected Output

Your final model should have:

- Complete atomic coordinates (N, CA, C, O, CB, side chains)
- Preserved secondary structure from scaffold
- Good stereochemistry (>95% Ramachandran favored)
- Low clash score (<5)
- Reasonable DOPE score (negative value)

## Summary

Congratulations! You've completed a full computational protein design pipeline:

1. **Generated** a novel protein backbone with RFdiffusion
2. **Designed** an amino acid sequence with ProteinMPNN
3. **Built** a complete atomic model with MODELLER

Your final model (`best_model.pdb`) represents a computationally designed two-helix bundle protein
ready for experimental validation or further computational analysis.

## Additional Resources

- MODELLER Documentation: https://salilab.org/modeller/
- PyMOL Tutorial: https://pymolwiki.org/
- Protein Structure Validation: https://saves.mbi.ucla.edu/
