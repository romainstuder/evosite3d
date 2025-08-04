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

### Step 1: Prepare Input Files

First, we need to prepare our scaffold and sequence:

```python
# prepare_inputs.py
from Bio import SeqIO, PDB
import os

def extract_ca_trace(pdb_file, output_file):
    """Extract only CA atoms from structure"""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('scaffold', pdb_file)

    class CASelect(PDB.Select):
        def accept_atom(self, atom):
            return atom.get_name() == 'CA'

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_file, CASelect())
    print(f"CA trace saved to {output_file}")

def prepare_alignment_file(sequence_file, template_pdb, output_ali):
    """Create alignment file for MODELLER"""
    # Read the designed sequence
    seq_record = list(SeqIO.parse(sequence_file, "fasta"))[0]
    sequence = str(seq_record.seq)

    # Get template sequence from PDB
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('template', template_pdb)

    template_seq = ""
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[0] == ' ':  # Standard residue
                    template_seq += 'X'  # Use X for unknown in template

    # Create alignment file
    with open(output_ali, 'w') as f:
        f.write(">P1;template\n")
        f.write("structure:template:FIRST:A:LAST:A::::\n")
        f.write(f"{template_seq}*\n\n")

        f.write(">P1;target\n")
        f.write("sequence:target::::::::\n")
        f.write(f"{sequence}*\n")

    print(f"Alignment file created: {output_ali}")

# Prepare files
extract_ca_trace('scaffold_clean.pdb', 'template_ca.pdb')
prepare_alignment_file('designed_sequence.fasta', 'template_ca.pdb', 'alignment.ali')
```

### Step 2: Create MODELLER Script

```python

```
