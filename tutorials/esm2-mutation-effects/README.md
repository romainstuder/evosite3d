# Predicting Mutation Effects on Protein Function with ESM-2

A hands-on tutorial for computational biologists.

[![Open Part 1 in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/romainstuder/evosite3d/blob/esm2-tutorial/tutorials/esm2-mutation-effects/01_zero_shot_scoring.ipynb) **Part 1: Zero-shot scoring**

---

## Why this tutorial?

Proteins evolve under competing pressures: they must remain stably folded while also performing
their function.
Ref:
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000002
https://pmc.ncbi.nlm.nih.gov/articles/PMC1360657/

In a previous study of RuBisCO, the enzyme that fixes CO₂ during photosynthesis,
illustrates this tension elegantly, showed that the evolution of C4 photosynthesis required
destabilizing mutations near the active site, which were tolerated only because earlier
stabilizing mutations had built up a structural "buffer," and were followed by compensatory
mutations that restored global stability. This stability-activity tradeoff shapes protein
evolution broadly, from natural adaptation to laboratory directed evolution. (Reference Studer
et al. (2014) https://www.pnas.org/doi/10.1073/pnas.1310811111)

In this tutorial, we explore the same phenomenon in a system with much richer experimental data:
**TEM-1 β-lactamase**, the enzyme responsible for bacterial resistance to penicillin-class
antibiotics. Firnberg et al. (2014) measured the fitness effect of nearly every possible single
mutation in TEM-1 by growing mutant libraries under ampicillin selection — giving us
experimental ground truth for ~5,000 variants, including classic examples of the tradeoff (G238S
extends the substrate range but destabilizes; M182T compensates).

**Our central question:** Can ESM-2, a protein language model trained only on sequences, predict
which mutations preserve, break, or enhance TEM-1 function?

## What you'll learn

By the end of this tutorial, you'll have:

- A working pipeline for mutation effect prediction with ESM-2
- An honest sense of where ESM-2 helps and where it doesn't
- A template you can adapt to your own protein of interest

## Prerequisites

- Familiarity with Python and PyTorch
- Basic understanding of transformers (helpful but not required)
- A laptop with 8GB+ RAM, **or** a free Google Colab account

Total runtime: ~1 hour for both notebooks combined.

## Tutorial structure

### Part 1 — Zero-shot scoring with ESM-2

uv run jupyter lab

**Notebook:** `01_zero_shot_scoring.ipynb`

We start without any training. ESM-2 was pretrained on ~65 million protein sequences using
masked language modeling, learning a rich prior over plausible protein sequences. We can exploit
this directly to score mutations.
https://en.wikipedia.org/wiki/Zero-shot_learning

You'll learn to:

- Load ESM-2 and explore its outputs
- Score mutations using **masked marginal likelihood**
- Compute a **per-position constraint profile** that reveals the active site without supervision
- Compare ESM-2 predictions to experimental fitness scores from Firnberg et al.

This part runs comfortably on any laptop (even an 8GB MacBook Air) in ~10 minutes.

**Key takeaway:** ESM-2 already "knows" a surprising amount about protein function from
pretraining alone. But its scores are relative, not calibrated to specific experimental conditions.

## How to run

### Option A: Google Colab (recommended for GPU access)

Click the "Open in Colab" badges at the top. The notebooks install dependencies automatically.
For Part 2, enable GPU via **Runtime → Change runtime type → T4 GPU**.

### Option B: Local installation

```bash
git clone https://github.com/romainstuder/evosite3d.git
cd evosite3d/tutorials/esm2-mutation-effects
uv sync
uv run jupyter lab
```

Then open the notebooks in order.

## Data

Both notebooks download the Firnberg et al. (2014) deep mutational scanning data from ProteinGym
automatically. No manual data preparation required.

## Going further

Once you've completed this tutorial, natural next steps include:

- **Frozen embeddings + classical ML** — extract ESM-2 embeddings once, train a ridge regression
  or gradient boosting model on top. Often within a few percent of full fine-tuning at a
  fraction of the compute.
- **Cross-protein generalization** — fine-tune on broader stability datasets (ProTherm, FireProt,
  MEGAscale) for models that transfer to new proteins.
- **Structure-aware models** — explore SaProt, ESM-3, or AlphaFold-derived embeddings when 3D
  structure information matters.
- **Multi-task learning** — predict stability, activity, and expression jointly from the same model.
- **Inverse folding** — use ESM-IF to design sequences for a given backbone structure.

## References

- **ESM-2:** Lin et al. (2023). _Evolutionary-scale prediction of atomic-level protein structure._ Science 379(6637):1123-1130.
- **TEM-1 DMS data:** Firnberg et al. (2014). _A comprehensive, high-resolution map of a gene's fitness landscape._ Mol Biol Evol 31(6):1581-1592.
- **RuBisCO motivation:** Studer et al. (2014). _Stability-activity tradeoffs constrain the adaptive evolution of RubisCO._ PNAS 111(6):2223-2228.
- **ProteinGym benchmark:** Notin et al. (2023). _ProteinGym: large-scale benchmarks for protein fitness prediction._ NeurIPS Datasets and Benchmarks.

## License

MIT

## Citation

If this tutorial helps your research or teaching, please cite or link back. Issues and pull
requests welcome.
