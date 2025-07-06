# EvoSite3D

A computational tool for analyzing evolutionary sites in 3D protein structures.

## Overview

EvoSite3D is designed to identify and analyze evolutionarily significant sites within protein structures by integrating sequence evolution data with 3D structural information. This tool helps researchers understand the relationship between protein evolution and structural constraints.

## Features

- **3D Structure Analysis**: Visualize evolutionary sites in protein 3D structures
- **Sequence Evolution Integration**: Combine phylogenetic analysis with structural data
- **Conservation Scoring**: Calculate conservation scores for residue positions
- **Structural Context**: Analyze evolutionary pressures in different structural environments
- **Visualization Tools**: Generate interactive 3D visualizations of evolutionary hotspots

## Installation

For detailed installation instructions, please refer to [INSTALL.md](INSTALL.md).

## Usage

### Basic Usage

```python
from evosite3d import EvoSite3D

# Initialize the analyzer
analyzer = EvoSite3D()

# Load protein structure and sequence alignment
analyzer.load_structure("protein.pdb")
analyzer.load_alignment("sequences.fasta")

# Analyze evolutionary sites
results = analyzer.analyze_sites()

# Generate 3D visualization
analyzer.visualize_sites(results, output="evolutionary_sites.png")
```

### Command Line Interface

```bash
# Analyze a single protein
evosite3d analyze --structure protein.pdb --alignment sequences.fasta --output results/

# Batch analysis
evosite3d batch --input_dir structures/ --alignments alignments/ --output results/

# Generate visualization
evosite3d visualize --results results/analysis.json --output visualization.png
```

## Input Data

### Required Files

1. **Protein Structure**: PDB format files containing 3D coordinates
2. **Sequence Alignment**: FASTA format multiple sequence alignment
3. **Phylogenetic Tree** (optional): Newick format tree file

### Supported Formats

- PDB, mmCIF for protein structures
- FASTA, CLUSTAL, PHYLIP for sequence alignments
- Newick, NEXUS for phylogenetic trees

## Output

The tool generates several types of output:

- **Conservation Scores**: Per-residue conservation analysis
- **Structural Context**: Analysis of evolutionary sites in different structural environments
- **3D Visualizations**: Interactive PyMOL sessions and static images
- **Statistical Reports**: Summary statistics and significance tests

## Examples

### Example 1: Single Protein Analysis

```python
# Analyze evolutionary sites in a kinase
analyzer = EvoSite3D()
analyzer.load_structure("kinase.pdb")
analyzer.load_alignment("kinase_family.fasta")

results = analyzer.analyze_sites(
    method="rate4site",
    structural_context=True,
    significance_threshold=0.05
)

# Identify catalytic sites under selection
catalytic_sites = results.filter_by_function("catalytic")
print(f"Found {len(catalytic_sites)} catalytic sites under selection")
```

### Example 2: Comparative Analysis

```python
# Compare evolution across protein families
families = ["kinases", "phosphatases", "proteases"]

for family in families:
    analyzer = EvoSite3D()
    analyzer.load_structure(f"{family}.pdb")
    analyzer.load_alignment(f"{family}_alignment.fasta")
    
    results = analyzer.analyze_sites()
    analyzer.save_results(f"{family}_evolution.json")
```

## Documentation

Detailed documentation is available in the `docs/` directory:

- [Installation Guide](docs/installation.md)
- [User Manual](docs/user_manual.md)
- [API Reference](docs/api_reference.md)
- [Tutorial](docs/tutorial.md)

## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Development Setup

```bash
# Clone with development dependencies
git clone https://github.com/romainstuder/evosite3d.git
cd evosite3d

# Install development dependencies
pip install -r requirements-dev.txt

# Run tests
pytest tests/

# Run linting
flake8 evosite3d/
black evosite3d/
```

## Citation

If you use EvoSite3D in your research, please cite:

```bibtex
@software{studer2024evosite3d,
  author = {Romain Studer},
  title = {EvoSite3D: Analyzing Evolutionary Sites in 3D Protein Structures},
  year = {2024},
  url = {https://github.com/romainstuder/evosite3d}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**Romain Studer**
- Computational Biologist at BenevolentAI
- Focus: Protein and nucleotide analysis
- LinkedIn: [romainstuder](https://linkedin.com/in/romainstuder)
- GitHub: [romainstuder](https://github.com/romainstuder)

## Acknowledgments

- Thanks to the computational biology community for feedback and testing
- Built using BioPython, PyMOL, and other excellent open-source tools
- Special thanks to [list any collaborators or institutions]

## Support

If you encounter any issues or have questions:

1. Check the [documentation](docs/)
2. Search existing [issues](https://github.com/romainstuder/evosite3d/issues)
3. Create a new issue with a minimal reproducible example

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history and updates.

---

*This project is under active development. Please report any bugs or feature requests through the GitHub issue tracker.*