# Installation Instructions

## Install Homebrew

On macOS, you can use homebrew: <https://brew.sh/> On Linux, you can use Linuxbrew:
https://docs.brew.sh/Homebrew-on-Linux

Or you could use Anaconda.

## Installation of common packages

```shell
brew install python \
             pyenv-virtualenv \
             uv
```

## Installation of bioinformatics packages

We will use a dedicated version of Homebrew for bioinformatics packages:
https://github.com/brewsci/homebrew-bio

```shell
brew tap brewsci/bio
brew install clustal-omega \
             dssp \
             fasttree \
             figtree \
             jalview \
             mafft \
             newick-utils \
             paml \
             phyml \
             pymol \
             trimal
```

### Links:

#### Package for ancestral sequence reconstruction (and many other things)

- CodeML / PAML: http://abacus.gene.ucl.ac.uk/software/paml.html

#### Alignment tools

- MAFFT: http://mafft.cbrc.jp/alignment/software/
- Clustal-Omega: http://www.clustal.org/omega/
- Clustal-Omega is the new aligner from ClustalW team, but much faster and more accurate
- TrimAl: <http://trimal.cgenomics.org/>

#### Alignment generation and visualisation

- Jalview: http://www.jalview.org/

#### Phylogenetics tools

- PhyML: http://www.atgc-montpellier.fr/phyml/binaries.php
- FastTree: http://www.microbesonline.org/fasttree/
- Newick utilities: <http://cegg.unige.ch/newick_utils>

# Tree visualisation

- FigTree: https://github.com/rambaut/figtree
- NJplot: http://doua.prabi.fr/software/njplot Note: I used NJplot a lot in the past. It is quick
  and simple, but it is now not compatible with recent macOS releases.

### Protein structure

- FoldX: <https://foldxsuite.crg.eu/>
- PyMOL: <https://www.pymol.org/>

## Installation of python modules

Python module dependencies (only `biopython` at the moment) are listed in the `pyproject.toml` file.
Just run to install:

```shell
uv pip install -r pyproject.toml
```

### Development Setup

```bash
# Clone with development dependencies
git clone https://github.com/romainstuder/evosite3d.git
cd evosite3d

# Create environments
uv venv
source .venv/bin/activate

# Install packages
uv pip install -r pyproject.toml

# (Optional) Install development dependencies
uv pip install -r pyproject.toml --extra dev
```

### Links:

- Biopython: http://biopython.org/wiki/Main_Page

## Scripts

The python scripts are in [./scripts/](./scripts/).

Make them accessible (i.e. in the $PATH list).

```shell
export PATH="$PATH:$HOME/Github/evosite3d/scripts";
chmod +x ./scripts/*.py
```
