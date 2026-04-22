# Positive Selection Detector

Detect pervasive positive selection on a protein-coding gene using
PAML CodeML site models (M8 vs M8a).

Given a gene symbol the pipeline:

1. Fetches the Ensembl Compara CDS alignment and gene tree.
2. Prunes to a curated panel of 41 reference vertebrates (fishes to
   mammals) so results are comparable across genes.
3. Runs CodeML **M8** and **M8a** in parallel and reports the
   likelihood-ratio test.
4. Extracts BEB positively-selected sites, maps them back to the
   reference protein and produces:
   - a per-site dN/dS plot (`_beb.png`),
   - a Jalview annotation file with BEB Pr(w>1) and dN/dS tracks
     (`_beb.jlv`),
   - a PyMOL script colouring significant sites on an AlphaFold or
     PDB structure (`_sites.pml`).

## Requirements

- Python >= 3.11
- [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) (`codeml` on PATH)
- [Newick Utilities](https://github.com/tjunier/newick_utils) (`nw_prune` on PATH)
- [trimAl](http://trimal.cgenomics.org/) (`trimal` on PATH, only needed
  when inputs are fetched from Ensembl)
- Python packages: `scipy`, `matplotlib`, `pandas`, `seaborn`
- Network access to `rest.ensembl.org` and `alphafold.ebi.ac.uk`

## Usage

```bash
# Basic run — resolves HLA-DQB1 via Ensembl, runs M8/M8a, full post-analysis
python run_codeml_site_models.py --gene-symbol HLA-DQB1

# Override UniProt accession and output directory
python run_codeml_site_models.py --gene-symbol HLA-DQB1 \
    --outdir HLA_DQB1_workdir --uniprot P01920

# Use the full vertebrate tree (fishes to mammals)
python run_codeml_site_models.py --gene-symbol TP53 --taxon 7742

# Use an experimental PDB structure instead of AlphaFold
python run_codeml_site_models.py --gene-symbol HLA-DQB1 \
    --pdb 1UVQ --chain A --resi-offset 32

# Only generate control files without running CodeML
python run_codeml_site_models.py --gene-symbol HLA-DQB1 --skip-codeml

# Stop after the LRT (no BEB extraction, plot or structure)
python run_codeml_site_models.py --gene-symbol HLA-DQB1 --skip-post-analysis
```

## Outputs

| File                     | Description                                        |
| ------------------------ | -------------------------------------------------- |
| `{prefix}_M8.ctl`        | CodeML control file for M8                         |
| `{prefix}_M8a.ctl`       | CodeML control file for M8a (null model)           |
| `{prefix}_M8.mlc`        | CodeML output for M8                               |
| `{prefix}_M8a.mlc`       | CodeML output for M8a                              |
| `{prefix}_M8.rst.txt`    | Preserved rst file from the M8 run                 |
| `{prefix}_beb_sites.tsv` | BEB sites with protein positions                   |
| `{prefix}_beb.txt`       | Cleaned per-site BEB probability table             |
| `{prefix}_beb.png`       | Per-site dN/dS Manhattan plot                      |
| `{prefix}_beb.jlv`       | Jalview annotation (BEB + dN/dS tracks)            |
| `{prefix}_sites.pml`     | PyMOL script highlighting BEB-significant residues |

## Reference species panel

The alignment is pruned to 41 vertebrate species spanning ~450 My of
divergence, grouped as: Primates (8), Rodentia (4), Lagomorpha (1),
Carnivora (2), Cetartiodactyla (5), Perissodactyla (1), Chiroptera (1),
Eulipotyphla (2), Marsupialia (1), Monotremata (1), Aves (3),
Reptilia (2), Amphibia (1), Fishes (9).
