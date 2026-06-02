# EvoSite3D

A comprehensive repository of scripts, software, and tutorials for molecular evolution, computational protein design, and structural bioinformatics.

## Overview

EvoSite3D integrates evolutionary sequence analysis, structural bioinformatics, and deep learning approaches to understand protein function, evolution, and design. This repository contains:

- **Production-ready tools** for evolutionary analysis and drug target validation
- **In-depth tutorials** covering sequence analysis, protein structure prediction and design, molecular dynamics, and computational biology
- **Research applications** analyzing viral evolution and genome datasets

Whether you're detecting positive selection, designing novel proteins, predicting mutation effects, or validating drug targets—this repository provides code, workflows, and educational materials to support your research.

## Quick Start

For detailed installation and setup, see [INSTALL.md](INSTALL.md).

All Python scripts and tools can be run using `uv run` (see [CLAUDE.md](CLAUDE.md) for environment details).

---

## 🛠️ Software Tools

### Target AI: AI-Powered Drug Target Validation

A Claude Code-powered assistant that validates therapeutic targets against diseases using real-time data from Open Targets Platform.

**Features:**

- Resolves gene symbols → Ensembl IDs and disease names → EFO IDs
- Fetches evidence from Open Targets (genetics, tractability, pathways, safety)
- Scores targets on Clinical Evidence, Druggability, and Pathway/Biology (0–5 scale each)
- Generates interactive HTML matrices and Markdown reports

**Location:** [`./software/target_ai/`](./software/target_ai/)

**Getting Started:**

```bash
# Install Claude Code
npm install -g @anthropic-ai/claude-code

# Navigate to the tool
cd ./software/target_ai
claude
```

Claude Code reads the embedded `CLAUDE.md` file and becomes your target validation assistant. No additional API keys needed—Open Targets Platform is free and open.

---

## 📚 Tutorials

### Evolutionary Analysis

#### [Detecting Pervasive Positive Selection (PAML Site Models)](./tutorials/pervasive_selection_site_models/)

Learn to detect continuous positive selection using CodeML/PAML. Covers dN/dS ratios, site-specific models (M0, M1a, M2a, M7, M8), and includes a complete HLA_DQB1 example with sequence alignment, phylogenetic trees, and structural visualisation.

#### [Detecting Episodic Positive Selection (Branch-Site Models)](./tutorials/positive_selection/)

Understand branch-site models for detecting positive selection in specific lineages. Includes real datasets and tree interpretation.

#### [Ancestral Sequence Reconstruction](./tutorials/ancestral_sequences_reconstruction/)

Reconstruct ancestral sequences at internal nodes of phylogenetic trees using maximum-likelihood methods.

### Protein Structure & Stability

#### [Estimating Mutation Stability Effects with FoldX](./tutorials/stability_impact/)

Learn to compute free energy changes (ΔG, ΔΔG) for protein mutations. Covers both FoldX 3 and FoldX 4 approaches. Suitable for understanding stability–activity trade-offs and mutation effects.

#### [Assessing Protein Structure Quality](./tutorials/assess_protein_structure_quality/)

Tools and workflows for validating predicted and experimental protein structures.

### Protein Design

#### [Computational Protein Design: Structure → Sequence → Atoms](./tutorials/generative-protein-design-workshop/)

A comprehensive three-part tutorial on de novo protein design:

1. **Structure Generation with RFdiffusion** — Generate novel protein backbones using diffusion models
2. **Sequence Design with ProteinMPNN** — Design amino acid sequences for your generated structures
3. **Full Atomic Modelling with MODELLER** — Build complete 3D models with all-atom resolution

Learn the full pipeline from backbone generation through atomic-level detail.

#### [Improving Protein-Protein Affinity](./tutorials/protein_design/improve_protein_protein_affinity/)

Design mutations to enhance protein-protein binding affinity using computational methods.

#### [Protein Design with RFdiffusion](./tutorials/generation_protein_rfdiffusion/)

Introduction to structure-based protein design using RFdiffusion for backbone generation.

### Mutation & Function Prediction

#### [Predicting Mutation Effects with ESM-2](./tutorials/esm2-mutation-effects/)

Use pre-trained protein language models (Facebook's ESM-2) to predict the functional impact of amino acid mutations without experimental training data. Covers the stability–activity trade-off in protein engineering.

#### [Predicting Steroid Receptor Selectivity with PyTorch](./tutorials/nuclear_receptor_selectivity/)

Build a multi-task graph neural network to predict small-molecule selectivity across six steroid hormone nuclear receptors (ERα, ERβ, AR, PR, GR, MR). Learn multi-task learning and latent space chemistry.

### Molecular Dynamics & Simulation

#### [Molecular Dynamics with OpenMM](./tutorials/molecular_dynamics/openmm_tutorial/)

Run production-ready molecular dynamics simulations of Ubiquitin using OpenMM. Includes system setup, equilibration, and analysis workflows.

### Omics & Systems Biology

#### [Graph Queries with Neo4j](./tutorials/omics/graph_queries_with_neo4j/)

Leverage graph databases to query and analyse biological networks and relationships.

### Post-Translational Modifications

#### [Phosphosite Identification with OpenMS](./tutorials/identification_phosposites/)

Detect and localise phosphorylation sites in proteins using mass spectrometry data and OpenMS workflows. Learn phosphoproteomics from data processing through site assignment.

### Foundations

#### [Learning PyTorch](./tutorials/learn_pytorch/)

Introduction to PyTorch for machine learning in computational biology, with practical examples.

#### [FDR & Q-Values](./tutorials/fdr_qvalue.md)

Statistical foundations for multiple testing correction and false discovery rates.

#### [NCBI Taxonomy](./tutorials/ncbi_taxonomy.md)

Working with NCBI taxonomy for phylogenetic and evolutionary studies.

---

## 🔬 Research

### Ebola Virus Glycoprotein Interface Mutations

**Location:** [`./research/2026_ebola_outbreak/`](./research/2026_ebola_outbreak/)

Comparative sequence analysis of Ebola virus glycoprotein (GP) between Zaire and Bundibugyo strains. Identifies 8 sequence differences at 21 antibody-contacting positions in the 3CSY crystal structure. Reveals regional clustering at the GP1-GP2 junction and conserved structural elements (disulfide bonds). Candidates for further experimental investigation.

---

## Key Features

- **3D Structure Integration** — Visualise evolutionary sites in protein structures
- **Evolutionary Analysis** — dN/dS ratios, positive selection detection, ancestral reconstruction
- **Mutation Prediction** — Predict stability and functional effects of mutations
- **Protein Design** — De novo design, sequence optimisation, affinity improvement
- **Molecular Simulation** — Molecular dynamics with OpenMM
- **Systems Biology** — Multi-omics integration and network analysis
- **AI-Powered Tools** — Leverage Claude Code for interactive target validation and analysis
- **Production Software** — Specialized tools for common computational biology tasks

## Working with Claude Code

Several tools in this repository are designed to work with **Claude Code**, Anthropic's official CLI for Claude:

1. **Target AI** (`./software/target_ai/`) — Start Claude Code in the directory; it reads `CLAUDE.md` and becomes your target validation assistant
2. **Interactive workflows** — Use Claude Code to explore data, refine analyses, and iterate on research questions

To get started with Claude Code:

```bash
npm install -g @anthropic-ai/claude-code
claude --help
```

---

## Installation & Setup

For detailed installation instructions, environment setup, and dependency management (using `uv`), see [INSTALL.md](INSTALL.md).

All Python code in this repository follows the guidelines in [CLAUDE.md](CLAUDE.md):

- Use `uv run` to execute Python scripts and tools
- Dependencies are managed via `uv add` / `uv remove`
- Tests run with `uv run pytest`
- Code is linted with `uv run ruff`

---

## Citation

If you use EvoSite3D in your research, please cite:

```bibtex
@software{studer2024evosite3d,
  author = {Romain Studer},
  title = {EvoSite3D: Analysing Evolutionary Sites in 3D Protein Structures},
  year = {2024},
  url = {https://github.com/romainstuder/evosite3d}
}
```

---

## Author

**Romain A. Studer**

- Senior Bioinformatics Data Scientist
- Previously affiliated with BenevolentAI, EMBL-EBI, UCL, and UNIL
- Focus: Protein and nucleotide analysis, computational biology, machine learning
- LinkedIn: [romainstuder](https://linkedin.com/in/romainstuder)
- GitHub: [romainstuder](https://github.com/romainstuder)

---

## Acknowledgements

- Computational biology community for feedback and testing
- Built with BioPython, PAML, PyMOL, OpenMM, ESM-2, RFdiffusion, and other excellent open-source tools
- Claude Code for interactive analysis and AI-assisted research workflows

---

## Support & Contributing

If you encounter issues or have questions:

1. **Search existing [issues](https://github.com/romainstuder/evosite3d/issues)** for solutions
2. **Create a new issue** with a minimal reproducible example and relevant details
3. **Contributions welcome** — pull requests, tutorials, and feedback are appreciated

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

_This project is under active development. Please report bugs, suggest features, or contribute improvements via the GitHub issue tracker._
