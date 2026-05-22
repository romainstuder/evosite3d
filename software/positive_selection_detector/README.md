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

| Step | Script                      | Role                                          |
| ---- | --------------------------- | --------------------------------------------- |
| 1    | `prepare_data.py`           | Fetch CDS alignment + tree + 3D structure     |
| 2    | `run_codeml_site_models.py` | Run CodeML M8 and M8a; LRT and BEB extraction |
| 2b   | `run_hyphy_site_models.py`  | Run HyPhy FUBAR (default) and/or BUSTED       |
| 3    | `analyze_site_models.py`    | dN/dS plot, Jalview annotation, PyMOL script  |

`main.nf` orchestrates all four under Nextflow with the obvious
parallelism.

## Requirements

- Python >= 3.11
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
python run_codeml_site_models.py --gene-symbol HLA-DQB1

# Step 2-bis: run HyPhy FUBAR (default) on the same alignment
python run_hyphy_site_models.py --gene-symbol HLA-DQB1

# Step 3: plots, Jalview annotation, PyMOL script
python analyse_site_models.py --gene-symbol HLA-DQB1
```

### Common options

```bash
# Use a different taxon (default 9347 = Boreoeutheria; 7742 = Vertebrata)
python prepare_data.py --gene-symbol TP53 --taxon 7742

# Override UniProt accession
python prepare_data.py --gene-symbol HLA-DQB1 --uniprot P01920

# Use an experimental PDB instead of AlphaFold
python prepare_data.py     --gene-symbol HLA-DQB1 --pdb 1UVQ
python analyse_site_models.py --gene-symbol HLA-DQB1 --pdb 1UVQ \
    --chain B --resi-offset 32

# Run a single CodeML model (designed for parallel execution)
python run_codeml_site_models.py --gene-symbol HLA-DQB1 --model M8
python run_codeml_site_models.py --gene-symbol HLA-DQB1 --model M8a
python run_codeml_site_models.py --gene-symbol HLA-DQB1 --extract-only

# HyPhy: choose method
python run_hyphy_site_models.py --gene-symbol HLA-DQB1 --method FUBAR
python run_hyphy_site_models.py --gene-symbol HLA-DQB1 --method BUSTED
python run_hyphy_site_models.py --gene-symbol HLA-DQB1 --method BOTH
```

### Nextflow

```bash
# Single gene
nextflow run main.nf --gene_symbols "HLA-DQB1"

# Multiple genes
nextflow run main.nf --gene_symbols "HLA-DQB1,TRIM5,FOXP2"

# Per-gene overrides via CSV
nextflow run main.nf --input genes.csv

# Skip HyPhy
nextflow run main.nf --gene_symbols "HLA-DQB1" --hyphy_method ""
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
| `{prefix}_busted.json`     | Raw BUSTED result (JSON)                            |
| `{prefix}_busted.log`      | HyPhy stdout for the BUSTED run                     |

### Visualise the results

```bash
# Jalview: open the protein MSA with the BEB / dN/dS annotation track
jalview -open HLA_DQB1/HLA_DQB1_subset.aa.mafft.fasta \
        -annotations HLA_DQB1/HLA_DQB1_beb.jlv

# PyMOL: load the structure and colour BEB-positive sites
pymol HLA_DQB1/HLA_DQB1_sites.pml
```

## Reference species panel

The alignment is pruned to 41 vertebrate species spanning ~450 My of
divergence, grouped as: Primates (8), Rodentia (4), Lagomorpha (1),
Carnivora (2), Cetartiodactyla (5), Perissodactyla (1), Chiroptera (1),
Eulipotyphla (2), Marsupialia (1), Monotremata (1), Aves (3),
Reptilia (2), Amphibia (1), Fishes (9).
