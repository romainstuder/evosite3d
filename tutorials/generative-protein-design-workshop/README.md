# Computational Protein Design: From Structure Generation to Full Atomic Models

## Introduction

Welcome to this comprehensive tutorial series on computational protein design. These practicals will
guide you through a complete protein design pipeline, combining state-of-the-art deep learning
methods with established structural biology tools.

### Overview of the Tutorial Series

This tutorial series consists of three interconnected practicals that demonstrate a modern approach
to de novo protein design:

1. **[Structure Generation with RFdiffusion](./01-structure-generation/README.md)** - Generate novel
   protein backbone structures using diffusion models
2. **[Sequence Design with ProteinMPNN](./02-sequence-design/README.md)** - Design amino acid
   sequences that fold into your generated structures
3. **[Full Atomic Modeling with MODELLER](./03-atomic-modelling/README.md)** - Build complete
   protein models with all atoms

### Learning Objectives

By completing these tutorials, you will:

- Understand the principles behind diffusion-based protein structure generation
- Learn how to use neural networks for protein sequence design
- Use homology modeling techniques for building full atomic structures
- Gain practical experience with cutting-edge protein design tools
- Develop skills to create novel proteins from scratch

### Prerequisites

To successfully complete these tutorials, you should have:

- Basic understanding of protein structure (primary, secondary, tertiary)
- Familiarity with command-line interfaces
- Python programming experience (intermediate level)
- Access to a GPU-enabled system (for RFdiffusion and ProteinMPNN)
- Basic knowledge of molecular biology concepts

### Required Software and Setup

Before starting, ensure you have the following tools installed:

1. **RFdiffusion** (v1.1.0 or later)

   - Repository: https://github.com/RosettaCommons/RFdiffusion
   - Requires: PyTorch, CUDA-enabled GPU
   - Note: macOS users could follow the installation there:
     https://github.com/romainstuder/RFdiffusion/blob/update-readme/README-macos.md

2. **ProteinMPNN** (latest version)

   - Repository: https://github.com/dauparas/ProteinMPNN
   - Requires: PyTorch, NumPy

3. **MODELLER** (v10.4 or later)

   - Website: https://salilab.org/modeller/
   - Requires: Academic license (free)

4. **PyMOL** (for visualisation)

   - Website: https://pymol.org/

5. **Python environment** with:
   - numpy, pandas, biopython
   - matplotlib, seaborn (for analysis)

### Tutorial Workflow

The workflow follows a logical progression from structure to sequence to full atomic model:

```
RFdiffusion → ProteinMPNN → MODELLER
     ↓              ↓            ↓
  Backbone   →  Sequence  → Full Atomic
  Structure     Design       Model
```

### What You'll Create

In this tutorial series, you will design a small two-helix bundle protein:

- **Size**: Approximately 40-60 amino acids
- **Structure**: Two alpha helices connected by a loop
- **Purpose**: Learn fundamental protein design principles with a manageable system

### Important Notes

- Each tutorial builds upon the previous one, so complete them in order
- Save all output files as they will be needed in subsequent steps
- Computational times vary based on your hardware
- Always validate your results using appropriate metrics

### Getting Help

If you encounter issues:

1. Check the official documentation for each tool
2. Verify your input file formats
3. Ensure all dependencies are correctly installed
4. Consult the troubleshooting section at the end of each tutorial

Let's begin with [Tutorial 1](./01-structure-generation/README.md): Generating a protein backbone
structure with RFdiffusion!
