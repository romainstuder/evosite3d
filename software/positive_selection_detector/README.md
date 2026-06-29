# Positive Selection Detector

Detect pervasive positive selection on a protein-coding gene using
PAML CodeML site models (M8 vs M8a), with HyPhy FUBAR/BUSTED as a
parallel cross-check.

Given a gene symbol the pipeline:

1. Fetches the Ensembl Compara CDS alignment and gene tree, plus a 3D
   structure (AlphaFold or RCSB PDB).
2. Prunes to a curated panel of 41 reference vertebrates (fishes to
   mammals) so results are comparable across genes.
3. Runs CodeML **M8** and **M8a** in parallel and reports the
   likelihood-ratio test; in parallel runs HyPhy **FUBAR** (or
   **BUSTED**) on the same alignment.
4. Extracts BEB positively-selected sites, maps them back to the
   reference protein and produces:
   - a per-site dN/dS plot (`_beb.png`),
   - a Jalview annotation file with BEB Pr(w>1) and dN/dS tracks
     (`_beb.jlv`),
   - a PyMOL script colouring significant sites on an AlphaFold or
     PDB structure (`_sites.pml`) — yellow spheres for CodeML BEB
     ≥ 0.95, orange sticks for BEB ≥ 0.50, cyan sticks for HyPhy
     FUBAR Pr ≥ 0.90.

## Pipeline scripts

The pipeline is split into focused steps so each can run in parallel
under Nextflow:

| Step | Script                             | Role                                              |
| ---- | ---------------------------------- | ------------------------------------------------- |
| 1    | `prepare_data.py`                  | Fetch CDS alignment + tree + 3D structure         |
| 2    | `run_codeml_site_models.py`        | Run CodeML M8 and M8a; LRT and BEB extraction     |
| 2b   | `run_hyphy_site_models.py`         | Run HyPhy FUBAR (default) and/or BUSTED           |
| 2c   | `run_codeml_branch_site_models.py` | Branch-site model A (MA vs MA0) along the lineage |
| 3    | `analyse_site_models.py`           | dN/dS plot, Jalview annotation, PyMOL script      |

`main.nf` orchestrates all four under Nextflow with the obvious
parallelism.

## Requirements

- Python >= 3.12
- [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) (`codeml` on PATH)
- [HyPhy](http://www.hyphy.org/) >= 2.5 (`hyphy` on PATH)
- [Newick Utilities](https://github.com/tjunier/newick_utils) (`nw_prune` on PATH)
- [trimAl](http://trimal.cgenomics.org/) (`trimal` on PATH)
- Python packages: `scipy`, `matplotlib`, `pandas`, `seaborn`, `biopython`
- Network access to `rest.ensembl.org`, `alphafold.ebi.ac.uk`,
  and `files.rcsb.org`
- Optional: [Nextflow](https://www.nextflow.io/) to run the full
  pipeline (`nextflow` on PATH)

## Usage

### Standalone (sequential)

```bash
# Step 1: fetch alignment, tree, and structure
uv run python prepare_data.py --gene-symbol HLA-DQB1

# Step 2: run both CodeML models, LRT, BEB extraction
uv run python run_codeml_site_models.py --gene-symbol HLA-DQB1

# Step 2-bis: run HyPhy FUBAR (default) on the same alignment
uv run python run_hyphy_site_models.py --gene-symbol HLA-DQB1

# Step 3: plots, Jalview annotation, PyMOL script
uv run python analyse_site_models.py --gene-symbol HLA-DQB1
```

### Common options

```bash
# Use a different taxon (default 9347 = Boreoeutheria; 7742 = Vertebrata)
uv run python prepare_data.py --gene-symbol TP53 --taxon 7742

# Override UniProt accession
uv run python prepare_data.py --gene-symbol HLA-DQB1 --uniprot P01920

# Use an experimental PDB instead of AlphaFold
uv run python prepare_data.py     --gene-symbol HLA-DQB1 --pdb 1UVQ
uv run python analyse_site_models.py --gene-symbol HLA-DQB1 --pdb 1UVQ \
    --chain B --resi-offset 32

# Run a single CodeML model (designed for parallel execution)
uv run python run_codeml_site_models.py --gene-symbol HLA-DQB1 --model M8
uv run python run_codeml_site_models.py --gene-symbol HLA-DQB1 --model M8a
uv run python run_codeml_site_models.py --gene-symbol HLA-DQB1 --extract-only

# HyPhy: choose method
uv run python run_hyphy_site_models.py --gene-symbol HLA-DQB1 --method FUBAR
uv run python run_hyphy_site_models.py --gene-symbol HLA-DQB1 --method BUSTED
uv run python run_hyphy_site_models.py --gene-symbol HLA-DQB1 --method BOTH
```

### HyPhy only

To run HyPhy without CodeML, first prepare the alignment (step 1), then
run the HyPhy step on its own:

```bash
# Step 1 (prerequisite): fetch the alignment HyPhy reads
uv run python prepare_data.py --gene-symbol HLA-DQB1

# HyPhy only (FUBAR by default; use --method to choose)
uv run python run_hyphy_site_models.py --gene-symbol HLA-DQB1
```

If step 1 has already been run for this gene, skip it and run
`run_hyphy_site_models.py` directly.

### Branch-site model A (lineage to the target)

`run_codeml_branch_site_models.py` runs CodeML branch-site **model A**
(`model = 2`, `NSsites = 2`) on every branch of the path from the root to
the target protein's leaf, testing each as the foreground (`#1`) branch.
For each branch it runs the alternative (MA) and null (MA0, `fix_omega = 1`,
`omega = 1`) models and reports the LRT (2ΔlnL vs a 50:50 mixture of a point
mass at 0 and χ²₁); BEB sites are extracted for branches that test
significant.

```bash
# Preview the lineage branches (no codeml)
uv run python run_codeml_branch_site_models.py --gene-symbol HLA-DQB1 --list-branches

# Run model A on every branch from the root to the target
uv run python run_codeml_branch_site_models.py --gene-symbol HLA-DQB1

# Run a single branch by index (1 = nearest the root); good for parallelism
uv run python run_codeml_branch_site_models.py --gene-symbol HLA-DQB1 --branch 3
```

Outputs (per branch `bNN`, numbered from the root towards the target):

| File                              | Description                                        |
| --------------------------------- | -------------------------------------------------- |
| `{prefix}_bNN.tree`               | Pruned tree with branch `bNN` marked `#1`          |
| `{prefix}_bsMA_bNN.mlc`           | CodeML output for the alternative model A          |
| `{prefix}_bsMA0_bNN.mlc`          | CodeML output for the null model A0                |
| `{prefix}_bsMA_bNN_beb_sites.tsv` | BEB sites with protein positions (significant)     |
| `{prefix}_branch_site_lrt.tsv`    | Per-branch LRT summary (2ΔlnL, p-value, BEB count) |

### Nextflow

```bash
# Single gene
nextflow run main.nf --gene_symbols "HLA-DQB1"

# Multiple genes
nextflow run main.nf --gene_symbols "HLA-DQB1,TRIM5,FOXP2"

# Per-gene overrides via CSV
nextflow run main.nf --input genes.csv

# Skip HyPhy (HyPhy is off by default; just omit --hyphy_method)
nextflow run main.nf --gene_symbols "HLA-DQB1"

# HyPhy only: skip CodeML (and the downstream analyse step) with the --skip_codeml flag
nextflow run main.nf --gene_symbols "HLA-DQB1" --skip_codeml --hyphy_method FUBAR

# Also run branch-site model A along the target lineage (off by default; slow)
nextflow run main.nf --gene_symbols "HLA-DQB1" --branch_site true
```

## Outputs

CodeML / post-analysis (under each gene workdir):

| File                     | Description                                        |
| ------------------------ | -------------------------------------------------- |
| `{prefix}_M8.ctl`        | CodeML control file for M8                         |
| `{prefix}_M8a.ctl`       | CodeML control file for M8a (null model)           |
| `{prefix}_M8.mlc`        | CodeML output for M8                               |
| `{prefix}_M8a.mlc`       | CodeML output for M8a                              |
| `{prefix}_M8.rst.txt`    | Preserved rst file from the M8 run                 |
| `{prefix}_beb_sites.tsv` | BEB sites with target-protein positions            |
| `{prefix}_beb.txt`       | Cleaned per-site BEB probability table             |
| `{prefix}_beb.png`       | Per-site dN/dS Manhattan plot                      |
| `{prefix}_beb.jlv`       | Jalview annotation (BEB + dN/dS tracks)            |
| `{prefix}_sites.pml`     | PyMOL script highlighting BEB-significant residues |

HyPhy:

| File                       | Description                                         |
| -------------------------- | --------------------------------------------------- |
| `{prefix}_fubar.json`      | Raw FUBAR result (JSON)                             |
| `{prefix}_fubar.log`       | HyPhy stdout for the FUBAR run                      |
| `{prefix}_fubar_sites.tsv` | FUBAR positive sites (Pr ≥ 0.90), with protein pos. |
| `{prefix}_fubar.jlv`       | Jalview annotation (FUBAR Pr + beta/alpha tracks)   |
| `{prefix}_busted.json`     | Raw BUSTED result (JSON)                            |
| `{prefix}_busted.log`      | HyPhy stdout for the BUSTED run                     |

`analyse_site_models.py` writes `_fubar.jlv` and adds FUBAR sites (cyan
sticks) to the PyMOL script. CodeML and HyPhy are handled independently,
so a HyPhy-only run still produces `_fubar.jlv` and a FUBAR-only `.pml`.

### Visualise the results

```bash
# Jalview: open the protein MSA with the CodeML BEB / dN/dS annotation track
jalview -open HLA_DQB1/HLA_DQB1_subset.aa.mafft.fasta \
        -annotations HLA_DQB1/HLA_DQB1_beb.jlv

# Jalview with the HyPhy FUBAR annotation track (HyPhy-only runs)
jalview -open HLA_DQB1/HLA_DQB1_subset.aa.mafft.fasta \
        -annotations HLA_DQB1/HLA_DQB1_fubar.jlv

# PyMOL: load the structure and colour positively-selected sites
pymol HLA_DQB1/HLA_DQB1_sites.pml
```

## Reference species panel

The alignment is pruned to 41 vertebrate species spanning ~450 My of
divergence, grouped as: Primates (8), Rodentia (4), Lagomorpha (1),
Carnivora (2), Cetartiodactyla (5), Perissodactyla (1), Chiroptera (1),
Eulipotyphla (2), Marsupialia (1), Monotremata (1), Aves (3),
Reptilia (2), Amphibia (1), Fishes (9).
