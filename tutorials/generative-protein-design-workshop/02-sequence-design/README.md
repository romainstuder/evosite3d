# Tutorial 2: Designing Amino Acid Sequences with ProteinMPNN

## Introduction

ProteinMPNN is a message-passing neural network that designs amino acid sequences for given protein
backbone structures. It learns the relationship between structure and sequence from known protein
structures, enabling it to propose sequences likely to fold into your target structure.

## Background: Inverse Folding with Neural Networks

### What is ProteinMPNN?

ProteinMPNN solves the "inverse folding" problem: given a 3D structure, what sequences will fold
into it? The model uses:

- **Graph Neural Networks**: Treats protein as a graph of interacting residues
- **Message Passing**: Propagates information between neighboring residues
- **Attention Mechanisms**: Focuses on structurally important interactions

### Key Advantages:

- Fast inference (~1 second per design)
- High success rate for designed sequences
- Handles multi-chain complexes
- Can incorporate sequence constraints

## Step-by-Step Tutorial

### Step 1: Prepare Input Structure

First, ensure your scaffold from Tutorial 1 is properly formatted:

```shell
export EVOSITE3D_SCRIPTS=./utils
python $EVOSITE3D_SCRIPTS/prepare_scaffold.py \
  --input_pdb=two_helix_scaffold.pdb \
  --output_pdb=scaffold_clean.pdb
```

### Step 2: Run ProteinMPNN

```shell
export MPNN_SCRIPTS=$HOME/Github/ProteinMPNN
```

```bash
# Basic sequence design
python $MPNN_SCRIPTS/protein_mpnn_run.py \
    --pdb_path scaffold_clean.pdb \
    --out_folder ./output/sequences \
    --num_seq_per_target 100 \
    --sampling_temp 0.1 \
    --batch_size 1
```

### Step 3: Understanding ProteinMPNN Output

The output includes:

- `sequences.fasta`: Designed sequences
- `scores.csv`: Confidence scores for each design

Example output format:

```
>scaffold_clean_A, score=0.8542, global_score=0.8542, T=0.1
MAEELKKLHEAVKLFEDVVLKHLKELVKLHEAAKLFEDVVL
>scaffold_clean_A, score=0.8321, global_score=0.8321, T=0.1
MSEELKKAHEAAKLFKDVVLKHLKKLVELHKAAKLFKDVVL
```

### Step 4: Analyse Sequence Properties

```shell
python $EVOSITE3D_SCRIPTS/analyse_sequences.py \
   --fasta_file=output/sequences/seqs/scaffold_clean.fa

```

### Step 6: Filter and Select Sequences

```shell
python $EVOSITE3D_SCRIPTS/select_sequences.py \
   --input_file=output/sequences/seqs/scaffold_clean.fa
```

### Step 7: Advanced ProteinMPNN Options

#### Design with Fixed Positions

```python
# Create position-specific design file
fixed_positions = {
    "scaffold_clean": {
        "A": [1, 5, 10, 15, 20, 25, 30, 35, 40]  # Fix these positions
    }
}

import json
with open('fixed_positions.json', 'w') as f:
    json.dump(fixed_positions, f)
```

```bash
# Run with constraints
python $MPNN_SCRIPTS/protein_mpnn_run.py \
    --pdb_path scaffold_clean.pdb \
    --fixed_positions_jsonl fixed_positions.json \
    --out_folder ./output/sequences_constrained \
    --num_seq_per_target 50
```

#### Temperature Sampling

```bash
# Lower temperature = more conservative designs
--sampling_temp 0.05  # Very conservative

# Higher temperature = more diverse designs
--sampling_temp 0.3  # More diverse

# Multiple temperatures
for temp in 0.05 0.1 0.2 0.3; do
    python run_ProteinMPNN.py \
        --pdb_path scaffold_clean.pdb \
        --out_folder ./output/sequences_T${temp} \
        --num_seq_per_target 20 \
        --sampling_temp ${temp}
done
```

### Step 8: Validate Sequence-Structure Compatibility

(You will need `mkdssp`. To install it on macOS:

```shell
brew tap brewsci/bio
brew install dssp
```

```shell
python $EVOSITE3D_SCRIPTS/validate_design.py \
   --sequences-file=output/sequences/seqs/scaffold_clean.fa \
   --structure-file=scaffold_clean.pdb
```

## Understanding ProteinMPNN Parameters

### Key Parameters:

1. **sampling_temp**: Controls sequence diversity

   - 0.05-0.1: Conservative, high confidence
   - 0.2-0.3: More diverse, exploratory
   - 1.0: Maximum diversity (rarely used)

2. **num_seq_per_target**: Number of sequences to generate

   - More sequences = better chance of success
   - Typical: 50-200 sequences

3. **model_name**: Different models available
   - `v_48_020`: Standard model
   - `v_48_030`: Soluble protein model
   - `v_48_010`: Monomeric protein model

## Common Issues and Solutions

### Issue 1: Low Diversity

```bash
# Increase temperature
--sampling_temp 0.2
# Generate more sequences
--num_seq_per_target 200
```

### Issue 2: Poor Scores

```bash
# Use soluble design model
--model_name v_48_030
# Check your structure quality
# python check_structure.py scaffold_clean.pdb
```

### Issue 3: Unusual Amino Acid Composition

```python
# Add composition constraints
bias_AA = {
    'A': [1.0] * len(sequence),  # Neutral bias
    'P': [0.0] * len(sequence),  # Avoid prolines
    'C': [0.1] * len(sequence),  # Limit cysteines
}
```

## Expected Output

Your designed sequences should have:

- High helix propensity (>1.0 average)
- Appropriate hydrophobic/hydrophilic balance
- Reasonable charge (-5 to +5)
- Good ProteinMPNN scores (>0.8)

Select your best sequence and save it as `designed_sequence.fasta` for Tutorial 3.

```shell
python $EVOSITE3D_SCRIPTS/extract_best_sequence.py \
    --input_fasta=filtered_sequences.fasta \
    --output_fasta=designed_sequence.fasta
cp scaffold_clean.pdb  ../03-atomic-modelling/scaffold_clean.pdb
cp designed_sequence.fasta  ../03-atomic-modelling/designed_sequence.fasta
```

## Next Steps

With your sequence designed, proceed to [Tutorial 3](../03-atomic-modelling/README.md) where you
will use MODELLER to build a full atomic model combining your scaffold structure and designed
sequence.

## Additional Resources

- ProteinMPNN Paper: [Dauparas et al., 2022](https://www.science.org/doi/10.1126/science.add2187)
- Official Repository: https://github.com/dauparas/ProteinMPNN
- Colab Notebook:
  https://colab.research.google.com/github/dauparas/ProteinMPNN/blob/main/colab_notebooks/quickstart.ipynb
