# Building Proteins Around Zinc-Finger Motifs with RFDiffusion

## Overview

Zinc-finger motifs are small structural domains that coordinate zinc ions and are commonly found in
DNA-binding proteins. This tutorial shows how to use RFDiffusion to design proteins that incorporate
zinc-finger motifs as structural scaffolds or functional elements.

## Zinc-Finger Motif Types

### Classical C2H2 Zinc Finger

- **Structure**: β-hairpin followed by α-helix
- **Coordination**: 2 cysteines, 2 histidines
- **Length**: ~30 residues
- **Sequence motif**: `X2-Cys-X2-4-Cys-X12-His-X3-5-His`

### C4 Zinc Finger (GATA-type)

- **Structure**: Two antiparallel β-sheets
- **Coordination**: 4 cysteines
- **Length**: ~25 residues

### RING Finger

- **Structure**: Cross-braced arrangement
- **Coordination**: Cys3His or Cys4His4
- **Length**: ~40-60 residues

## Method 1: Motif Scaffolding with Existing Zinc-Finger

### Step 1: Prepare Input Structure

First, extract or create a zinc-finger motif PDB file:

```shell
# Download a zinc-finger protein (example: 1ZAA - classic C2H2)
wget https://files.rcsb.org/download/1ZAA.pdb
```

```shell
# Optional: visualise it
pymol 1ZAA.pdb
```

```shell
# Extract just the zinc-finger domain (residues 11-39)
rfd_extract_motif.py --input 1ZAA.pdb --chain A --residues 11-39 --output zf_motif.pdb
```
```shell
# Optional: visualise it
pymol zf_motif.pdb
```


**Python script to extract motif (extract_motif.py):**


### Step 2: Design Scaffold Around Motif

```bash
pyenv local 3.10.16
python $RFD_SCRIPTS/run_inference.py \
    inference.output_prefix=zf_scaffold \
    inference.input_pdb=zf_motif.pdb \
    inference.num_designs=10 \
    inference.ckpt_override_path=./models/Base_ckpt.pt \
    '+contigs=[20-40/A11-39/20-40]' \
    '+hotspot_res=[A14,A17,A32,A35]'
```

**Parameters explained:**

- `20-40/A11-39/20-40`: Design 20-40 residues before motif, keep motif (A11-39), design 20-40 after
- `hotspot_res`: Zinc-coordinating residues (Cys14, Cys17, His32, His35)
- Generates scaffolds of 71-119 total residues

## Method 2: De Novo Design with Zinc-Finger Constraints

### Step 1: Create Zinc-Finger Template

```shell
rfd_create_zf_template.py
``

### Step 2: Design with Template Constraints

```bash
python run_inference.py \
    inference.output_prefix=denovo_zf \
    inference.input_pdb=zf_template.pdb \
    '+contig_map.contigs=[10-20/A1-4/15-25/A1-4/10-20]' \
    '+ppi.hotspot_res=[A1,A2,A3,A4]' \
    inference.num_designs=15 \
    inference.ckpt_override_path=./models/Base_ckpt.pt \
    potentials.guiding_potentials=["type:olig_contacts,weight_intra:1,weight_inter:0.1"] \
    potentials.guide_scale=2.0 \
    potentials.guide_decay="quadratic"
```

## Method 3: Multi-Finger Domain Design

Design a protein with multiple zinc fingers (like transcription factors):

### Step 1: Prepare Multi-Finger Template

```bash
# Create a 3-finger domain template
python create_multi_zf.py --fingers 3 --spacing 7 --output triple_zf.pdb
```

```python
# create_multi_zf.py
import argparse
import numpy as np
from Bio.PDB import Structure, Model, Chain, Residue, Atom, PDBIO

def create_multi_zf(num_fingers, finger_spacing, output_file):
    structure = Structure.Structure('Multi_ZF')
    model = Model.Model(0)
    chain = Chain.Chain('A')

    # Standard C2H2 zinc finger template (30 residues)
    zf_template = """
    PYKCELCGKSFRSDSSLTKHQRIHTGEKPFA
    """

    residue_id = 1
    x_offset = 0

    for finger in range(num_fingers):
        for i, aa in enumerate(zf_template.strip()):
            # Simple helical geometry for demonstration
            x = x_offset + i * 1.5
            y = 5 * np.sin(i * 0.3)
            z = 5 * np.cos(i * 0.3)

            residue = Residue.Residue((' ', residue_id, ' '), aa, ' ')
            atom = Atom.Atom('CA', [x, y, z], 0.0, 1.0, ' ', 'CA', residue_id)
            residue.add(atom)
            chain.add(residue)
            residue_id += 1

        x_offset += 45  # Space between fingers

        # Add linker residues
        if finger < num_fingers - 1:
            linker = "TGEKPFA"  # 7-residue linker
            for i, aa in enumerate(linker):
                x = x_offset + i * 1.5
                y = 0
                z = 0

                residue = Residue.Residue((' ', residue_id, ' '), aa, ' ')
                atom = Atom.Atom('CA', [x, y, z], 0.0, 1.0, ' ', 'CA', residue_id)
                residue.add(atom)
                chain.add(residue)
                residue_id += 1

    model.add(chain)
    structure.add(model)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fingers', type=int, default=3)
    parser.add_argument('--spacing', type=int, default=7)
    parser.add_argument('--output', default='multi_zf.pdb')
    args = parser.parse_args()

    create_multi_zf(args.fingers, args.spacing, args.output)
```

### Step 2: Design Multi-Finger Protein

```bash
python run_inference.py \
    inference.output_prefix=multi_zf_protein \
    inference.input_pdb=triple_zf.pdb \
    'contig_map.contigs=[20-40/A1-97/20-40]' \
    'ppi.hotspot_res=[A3,A6,A23,A26,A33,A36,A53,A56,A63,A66,A83,A86]' \
    inference.num_designs=20 \
    inference.ckpt_override_path=./models/Base_ckpt.pt \
    potentials.guiding_potentials=["type:olig_contacts,weight_intra:1"] \
    potentials.guide_scale=1.5
```

## Working Example: Single Zinc-Finger Scaffold

Let's create a complete working example:

### Step 1: Download and Prepare Zinc-Finger Structure

```bash
# Get a well-characterized zinc finger (Zif268)
wget https://files.rcsb.org/download/1ZAA.pdb

# Extract finger 1 (residues 32-61)
python -c "
from Bio.PDB import PDBParser, PDBIO, Select

class ZFSelect(Select):
    def accept_residue(self, residue):
        return 32 <= residue.id[1] <= 61 and residue.parent.id == 'A'

parser = PDBParser(QUIET=True)
structure = parser.get_structure('zif268', '1ZAA.pdb')
io = PDBIO()
io.set_structure(structure)
io.save('zf_finger1.pdb', ZFSelect())
print('Extracted zinc finger motif to zf_finger1.pdb')
"
```

### Step 2: Design Scaffold

```bash
python run_inference.py \
    inference.output_prefix=zf_scaffold_example \
    inference.input_pdb=zf_finger1.pdb \
    'contig_map.contigs=[15-30/A32-61/15-30]' \
    'ppi.hotspot_res=[A35,A38,A52,A55]' \
    inference.num_designs=5 \
    inference.ckpt_override_path=./models/Base_ckpt.pt \
    denoiser.noise_scale_ca=1.0 \
    denoiser.noise_scale_frame=1.0 \
    inference.ckpt_override_path=./models/Base_ckpt.pt
```

### Step 3: Analyze Results

```python
# analyze_zf_designs.py
import sys
from Bio.PDB import PDBParser
import numpy as np

def analyze_zf_design(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('design', pdb_file)

    # Find zinc-coordinating residues
    zn_coords = []
    for residue in structure.get_residues():
        if residue.resname in ['CYS', 'HIS']:
            if residue.resname == 'CYS' and 'SG' in residue:
                zn_coords.append(residue['SG'].coord)
            elif residue.resname == 'HIS' and 'NE2' in residue:
                zn_coords.append(residue['NE2'].coord)

    if len(zn_coords) >= 4:
        # Calculate potential zinc position (centroid)
        zn_pos = np.mean(zn_coords[:4], axis=0)

        # Calculate distances to coordinating atoms
        distances = [np.linalg.norm(coord - zn_pos) for coord in zn_coords[:4]]

        print(f"Design: {pdb_file}")
        print(f"Zinc coordination distances: {distances}")
        print(f"Average distance: {np.mean(distances):.2f} Å")
        print(f"Distance std: {np.std(distances):.2f} Å")
        print("Good zinc coordination: 2.0-2.5 Å average")
        print()

# Analyze all designs
for i in range(5):
    analyze_zf_design(f'zf_scaffold_example_{i}.pdb')
```

### Step 4: Validate Structure Quality

```bash
# Basic structure validation
python -c "
import biotite.structure.io as bsio
import biotite.structure as struc

for i in range(5):
    pdb_file = f'zf_scaffold_example_{i}.pdb'
    try:
        structure = bsio.load_structure(pdb_file)

        # Calculate radius of gyration
        ca_atoms = structure[structure.atom_name == 'CA']
        rg = struc.gyration_radius(ca_atoms)

        # Count secondary structure (simplified)
        length = len(ca_atoms)

        print(f'Design {i}: Length={length}, Rg={rg:.1f}Å')
    except:
        print(f'Could not analyze design {i}')
"
```

## Post-Processing and Optimization

### 1. Sequence Design with ProteinMPNN

```bash
# Install ProteinMPNN if not available
git clone https://github.com/dauparas/ProteinMPNN.git

# Design sequences for the best scaffold
python ProteinMPNN/protein_mpnn_run.py \
    --pdb_path zf_scaffold_example_0.pdb \
    --pdb_path_chains A \
    --fixed_positions "35,38,52,55" \
    --num_seq_per_target 10 \
    --out_folder sequences/
```

### 2. Structure Prediction Validation

```bash
# Use ColabFold for quick structure prediction
python -c "
import subprocess
import os

# Install colabfold if needed
subprocess.run(['pip', 'install', 'colabfold[alphafold]'])

# Predict structures for designed sequences
for seq_file in os.listdir('sequences/'):
    if seq_file.endswith('.fa'):
        cmd = f'colabfold_batch sequences/{seq_file} predictions/ --num-models 1'
        subprocess.run(cmd.split())
"
```

## Tips for Zinc-Finger Design Success

### 1. Preserve Coordination Geometry

- Keep zinc-coordinating residues (Cys, His) as hotspots
- Maintain proper spacing between coordinating residues
- Typical Zn-S distance: 2.3 Å, Zn-N distance: 2.1 Å

### 2. Consider Fold Stability

- Include sufficient secondary structure elements
- Ensure proper hydrophobic core formation
- Add stabilizing residues in scaffold regions

### 3. Functional Considerations

- For DNA-binding: preserve recognition helix orientation
- For protein-protein interactions: maintain binding interface
- Consider metal accessibility and binding kinetics

### 4. Validation Pipeline

1. **Structural validation**: Check coordination geometry
2. **Sequence design**: Use ProteinMPNN with constraints
3. **Fold prediction**: Validate with AlphaFold2/ColabFold
4. **Experimental testing**: Express, purify, test metal binding

## Troubleshooting Common Issues

**Problem**: Poor zinc coordination geometry **Solution**: Increase hotspot weights, reduce noise
scales

**Problem**: Unstable folds **Solution**: Add more scaffold residues, use stability potentials

**Problem**: Loss of function **Solution**: Include more functional residues as hotspots

**Problem**: Metal binding sites too buried **Solution**: Design with surface exposure constraints

This comprehensive approach allows you to design stable, functional proteins incorporating
zinc-finger motifs for various applications including transcription factors, biosensors, and
structural scaffolds.
