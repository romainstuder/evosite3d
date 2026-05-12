# Predicting Steroid Receptor Selectivity with PyTorch

A tutorial on training a multi-task graph neural network to predict whether a small molecule
binds to any of the six **steroid hormone nuclear receptors**:

- **ERα** (estrogen receptor alpha) — endogenous ligand: estradiol
- **ERβ** (estrogen receptor beta) — endogenous ligand: estradiol
- **AR** (androgen receptor) — endogenous ligand: testosterone / DHT
- **PR** (progesterone receptor) — endogenous ligand: progesterone
- **GR** (glucocorticoid receptor) — endogenous ligand: cortisol
- **MR** (mineralocorticoid receptor) — endogenous ligand: aldosterone

The goal is to learn what selectivity actually looks like in latent space across the full subfamily.

## Files

- `README.md` — this file (concepts and background)
- `nuclear_receptor_selectivity.ipynb` — runnable notebook with all code

## Why six receptors makes this richer

These six proteins form the NR3 subfamily and split into two groups:

- **NR3A — estrogen receptors:** ERα and ERβ. Phenolic A-ring ligands. Smaller, more hydrophobic pocket.
- **NR3C — 3-ketosteroid receptors:** AR, PR, GR, MR. Bind 4-en-3-one steroid scaffolds. Larger pocket.

Within each group there is real cross-reactivity, and across groups there is little. Concrete examples:

- ERα vs ERβ: pockets differ by only two residues; many ligands hit both. Selective ligands (DPN for ERβ, PPT for ERα are pharmacologically valuable.
- AR vs PR: progesterone is a partial AR ligand; many AR antagonists hit PR.
- GR vs MR: pockets are extremely similar; spironolactone (an MR antagonist) hits AR and PR; aldosterone is the only
- well-discriminated MR ligand.
- ER vs the 3-ketosteroid group: very little overlap.

A model that has actually learned chemistry should reproduce this block structure: ERα and ERβ as one cluster,
AR/PR/GR/MR as another, with cross-block confusion being rare.

## What you will build

A graph convolutional network with a shared trunk and **six** output heads — one per receptor.

```
SMILES ──► molecular graph ──► GCN trunk ──► shared MLP ──┬─► ERα head
                                                          ├─► ERβ head
                                                          ├─► AR head
                                                          ├─► PR head
                                                          ├─► GR head
                                                          └─► MR head
```

## Key concepts

### 1. Graph neural networks for molecules

Small molecules are naturally graphs. Each atom is a node carrying features (element,
hybridization, formal charge, aromaticity, hydrogen count). Each bond is an edge. A GCN
propagates information between neighboring atoms over several layers, then pools all atom
embeddings into a single graph-level vector. This is more expressive than fixed fingerprints
because the network learns which substructures matter for the task.

### 2. Multi-task learning with masked loss

Most molecules in ChEMBL are measured against only one or two receptors. ERα and AR are the most-tested; MR is the least. A naive approach would drop molecules with missing labels and lose almost all the data.

The fix is a **masked binary cross-entropy loss**. For each molecule we compute six logits and
only count the loss on receptors where we actually have a measurement. Missing labels are tagged with `-1` and zeroed out in the loss:

```python
mask = (y != -1).float()
loss = bce(logits, y.clamp(min=0)) * mask
return loss.sum() / mask.sum()
```

Six receptors makes the masked-loss benefit even more pronounced: the shared trunk sees every molecule regardless of which receptor it was tested on, and learns a single representation of "steroid-like chemistry" that all heads draw on.

### 3. What "binder" means

Bioactivity comes as IC50, Ki, Kd, or EC50 in nM. We convert to a pActivity (`-log10` of molar concentration) and call anything with pActivity ≥ 6 (≤ 1 µM) a binder. This is a coarse threshold but standard in cheminformatics, and the masked-loss machinery generalizes if you want to do regression instead.

### 4. The selectivity story

When you embed test molecules with the trained network and run t-SNE, expect to see:

- A clear **ER cluster** (ERα + ERβ ligands together) — phenolic, often non-steroidal scaffolds like raloxifene, tamoxifen, genistein
- A larger **3-ketosteroid cluster** (AR + PR + GR + MR ligands) — classic steroid scaffolds
- Within the 3-ketosteroid blob, sub-structure: AR ligands somewhat distinct (DHT-like), GR/MR strongly overlapping, PR sitting in between

The ERα/ERβ ligands will be hard to separate from each other within the ER cluster — that is the model honestly representing the biology, not a bug. Same for GR/MR.

## Expected results

After 50 epochs on a typical ChEMBL pull, AUROCs roughly follow data volume and pocket distinctiveness:

| Receptor | Expected AUROC | Notes                          |
| -------- | -------------- | ------------------------------ |
| ERα      | 0.88–0.93      | Most data, distinct pocket     |
| AR       | 0.85–0.90      | Lots of prostate cancer data   |
| PR       | 0.80–0.85      | Moderate data                  |
| GR       | 0.80–0.88      | Moderate data                  |
| ERβ      | 0.75–0.85      | Less data, very similar to ERα |
| MR       | 0.65–0.75      | Least data, very similar to GR |

A scaffold split (instead of random) typically gives a more conservative estimate of generalization,
since structurally related molecules are forced into different splits. AUROCs commonly drop by
roughly 0.05–0.10 under that setup, though the exact gap depends on the dataset and on the
scaffold-clustering choices. The tutorial uses a random split for simplicity.

## The cross-reactivity matrix

A nice diagnostic with six receptors: compute Pearson correlations between predicted probabilities across receptors. You should see:

- **High correlation within blocks** (ERα↔ERβ, GR↔MR, AR↔PR predictions correlate)
- **Low correlation across blocks** (ER vs 3-ketosteroid predictions decorrelate)

The notebook computes and plots this 6×6 matrix. The NR3A / NR3C block structure should pop out
clearly — the model reconstructs the subfamily split purely from ligand chemistry.

## Running on Apple Silicon

The notebook runs comfortably on an M1 MacBook with 16 GB RAM:

- Dataset: ~10–20k unique molecules across all six receptors
- Model: ~55k parameters
- Training: 5–10 minutes on CPU; MPS optional
- Peak RAM: well under 1 GB

One slow step is the initial ChEMBL API pull (10–15 minutes for six targets, network-bound).
The notebook caches each target separately to parquet so you can resume if it's interrupted.

## Setup

Instructions to install dependencies are [here](../../INSTALL.md)

```bash
cd evosite3d/tutorials/nuclear_receptor_selectivity
uv sync
uv run jupyter lab nuclear_receptor_selectivity.ipynb
```

If `torch-geometric` complains about `torch-scatter` or `torch-sparse`, you usually don't need
them for this tutorial — the basic `GCNConv` works without them.
