# README

## Install Homebrew
On macOS, you can use homebrew: <https://brew.sh/>
On Linux, you can use Linuxbrew: https://docs.brew.sh/Homebrew-on-Linux

Or you could you use Anaconda.

## Installation of common packages

```shell
brew install python \
             pyenv-virtualenv \
             uv
```

Links:
- BioPython: <http://biopython.org/>
- Jalview: <http://www.jalview.org/>
- MAFFT: <http://mafft.cbrc.jp/alignment/software/>
- Newick utilities: <http://cegg.unige.ch/newick_utils>
- PAML: <http://abacus.gene.ucl.ac.uk/software/paml.html>
- PyMOL: <https://www.pymol.org/>
- R: <https://www.r-project.org/>
- TrimAl: <http://trimal.cgenomics.org/>


## Installation of bioinformatics packages

We will use a dedicated version of Homebrew for bioinformatics packages:
https://github.com/brewsci/homebrew-bio

```shell
brew tap brewsci/bio
brew install clustal-omega \
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

Two packages are currently problematic, and would need to be installed from sources:
- Newick-utils: https://github.com/tjunier/newick_utils
- PhyML: https://github.com/stephaneguindon/phyml/




Links:
```shell
# Package for ancestral sequence reconstruction (and many other things)
- CodeML from PAML: http://abacus.gene.ucl.ac.uk/software/paml.html

# Alignment tools
- MAFFT L-INS-i: http://mafft.cbrc.jp/alignment/software/
- Clustal-Omega: http://www.clustal.org/omega/
- Clustal-Omega is the new aligner from ClustalW team, but much faster and more accurate

# Alignment visualisation
- Jalview: http://www.jalview.org/

# Phylogenetics tools
- PhyML: http://www.atgc-montpellier.fr/phyml/binaries.php
- FastTree: http://www.microbesonline.org/fasttree/

```

Note: I used NJplot a lot  in the past. It is quick and simple, but it is now not compatible
with recent macOS releases.
# Tree visualisation
- NJplot: http://doua.prabi.fr/software/njplot


## Installation of python modules

Python module dependencies (only `biopython` at the moment) are listed in the `pyproject.toml`
file. Just run to install:
```shell
uv pip install -r pyproject.toml
```

Links:
- Biopython: http://biopython.org/wiki/Main_Page


## Scripts

The python scripts you need to install are in <https://github.com/romainstuder/evosite3d/scripts/>:


Make any accessible directory (i.e. in the working directory or in the $PATH list).

```shell
export PATH="$PATH:$HOME/Github/evosite3d/scripts";
chmod +x ./scripts/*.py
```
