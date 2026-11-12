# Installation Instructions

## Install Homebrew

On macOS, you can use homebrew: <https://brew.sh/>.
On Linux, you can use Linuxbrew: <https://docs.brew.sh/Homebrew-on-Linux>.

## Installation of common packages

```shell
brew install python uv
```

## Installation of bioinformatics packages

We will use a dedicated version of Homebrew for bioinformatics packages:
<https://github.com/brewsci/homebrew-bio>.

```shell
brew tap brewsci/bio
brew install clustal-omega \
             dssp \
             fasttree \
             hyphy \
             jalview \
             mafft \
             newick-utils \
             paml \
             phyml \
             pymol \
             trimal
```

(Note: figtree has been removed from brewsci/bio)

### Links:

#### Package for ancestral sequence reconstruction (and many other things)

- CodeML / PAML: <http://abacus.gene.ucl.ac.uk/software/paml.html>

#### Alignment tools

- MAFFT: <http://mafft.cbrc.jp/alignment/software/>
- Clustal-Omega: <http://www.clustal.org/omega/> — the successor to ClustalW from the same team,
  much faster and more accurate
- TrimAl: <http://trimal.cgenomics.org/>

#### Alignment generation and visualisation

- Jalview: <http://www.jalview.org/>

#### Phylogenetics tools

- PhyML: <http://www.atgc-montpellier.fr/phyml/binaries.php>
- FastTree: <http://www.microbesonline.org/fasttree/>
- Newick utilities: <http://cegg.unige.ch/newick_utils>

#### Tree visualisation

- FigTree: <https://github.com/rambaut/figtree>
- NJplot: <http://doua.prabi.fr/software/njplot> — I used NJplot a lot in the past; quick and
  simple, but no longer compatible with recent macOS releases.

#### Protein structure

- FoldX: <https://foldxsuite.crg.eu/>
- PyMOL: <https://www.pymol.org/>

## Installation of Python modules

Python dependencies are declared in `pyproject.toml` and locked in `uv.lock`. uv manages the
virtual environment automatically — do not create or activate `.venv` yourself.

```shell
# Pin the interpreter the project expects (one-time)
uv python install 3.12

# Install runtime dependencies into a managed .venv
uv sync
```

Run anything Python-related through `uv run`:

```shell
uv run jupyter lab
uv run pytest
uv run ruff check .
```

### Development setup

```shell
git clone https://github.com/romainstuder/evosite3d.git
cd evosite3d

# Install runtime + dev dependencies (PEP 735 dependency group)
uv sync --group dev

# Install pre-commit hooks (use uvx so it doesn't leak into the project venv)
uvx pre-commit install
```

### Python ecosystem links

- Biopython: <http://biopython.org/wiki/Main_Page>

## Specialised tools for individual tutorials

The dependencies above cover the core evolution, structure, machine-learning, and
molecular-dynamics tutorials. A few tutorials rely on heavier or external tools that are not
installed by `uv sync` — install these only if you intend to run the corresponding tutorial.

### Claude Code (software/target_ai)

The Target AI drug-target validation tool runs through Claude Code, Anthropic's official CLI:

```shell
npm install -g @anthropic-ai/claude-code
```

You will need an Anthropic API key or a Claude Pro/Max subscription. See
<https://docs.anthropic.com/en/docs/claude-code/overview>.

### OpenMS (tutorials/identification_phosposites)

Phosphosite identification uses the OpenMS mass-spectrometry suite. Install the command-line tools
and the Python bindings:

```shell
brew install openms      # command-line tools
uv pip install pyopenms  # Python bindings (into the managed .venv)
```

See <https://openms.de/>.

### Neo4j (tutorials/omics/graph_queries_with_neo4j)

Graph queries require a running Neo4j database. Install Neo4j Desktop or the community server from
<https://neo4j.com/download/>

### Protein design workshop (tutorials/generative-protein-design-workshop, tutorials/generation_protein_rfdiffusion)

The de novo design tutorials use GPU-based deep-learning tools with their own (conda/git) install
procedures. Follow the upstream instructions and the per-tutorial READMEs:

- RFdiffusion: <https://github.com/RosettaCommons/RFdiffusion>
- ProteinMPNN: <https://github.com/dauparas/ProteinMPNN>
- MODELLER: <https://salilab.org/modeller/> (academic licence required)

### FoldX (tutorials/stability_impact)

FoldX is not available through Homebrew. Download the binary after registering for a free academic
licence at <https://foldxsuite.crg.eu/>.

## Scripts

The Python scripts are in [./scripts/](./scripts/).

Make them accessible (i.e. add them to `$PATH`):

```shell
export PATH="$PATH:$HOME/Github/evosite3d/scripts"
chmod +x ./scripts/*.py
```
