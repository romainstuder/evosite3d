# Predicting Mutation Effects on Protein Function with ESM-2

A hands-on tutorial for computational biologists.

---

## Why this tutorial?

### The stability–activity tradeoff

Proteins live under two competing pressures. They must remain **stably folded** in their cellular
environment, and they must **perform a function** — bind a substrate, catalyze a reaction, signal
to a partner. These two requirements pull in opposite directions at the residue level:

- **Stability** is a global property. It comes from a network of well-packed hydrophobic cores,
  satisfied hydrogen bonds, and favorable electrostatics distributed across the fold. Most random
  mutations chip away at this network — the average single mutation costs roughly +1 kcal/mol of
  unfolding free energy ([Tokuriki et al. 2008][tokuriki]).
- **Activity** is local and chemically demanding. Catalytic residues are typically charged or
  polar groups placed in geometrically strained, often desolvated positions — exactly the
  configurations the folding energy "wants" to avoid. The same residues that make the active site
  reactive also make it a stability liability ([Chen et al. 2005][chen] dissect this for the
  CTX-M β-lactamase, where catalytic Lys73 and Glu166 sit in a strained, water-mediated
  geometry).

The consequence is a **stability–activity tradeoff**: mutations that introduce or refine activity
are frequently _destabilizing_, and mutations that stabilize the fold are often _neutral or
slightly deleterious for activity_. Tokuriki and Tawfik analyzed 548 function-altering mutations
across 22 enzymes and found them, on average, to be no more destabilizing than random mutations —
roughly +0.9 kcal/mol — but this is enough that, without compensation, the protein would unfold.
What lets evolution navigate this constraint is the accumulation of apparently _silent_
**compensatory mutations** that buy back stability ahead of, alongside, or after the
function-altering change.

During my postdoc I worked on RuBisCO — the enzyme that fixes CO₂ in photosynthesis — and
observed the same pattern there. The evolution of C4 photosynthesis required destabilizing
mutations near the active site, which were tolerated only because earlier stabilizing mutations
had built up a structural "buffer," and were followed by compensatory mutations that restored
global stability ([Studer et al. 2014][studer]). This pattern shapes protein evolution broadly,
from natural adaptation to laboratory directed evolution.

### TEM-1 β-lactamase as a model system

In this tutorial, we explore the same phenomenon in a system with much richer experimental data:
**TEM-1 β-lactamase**, the enzyme responsible for bacterial resistance to penicillin-class
antibiotics. Firnberg et al. (2014) measured the fitness effect of nearly every possible single
mutation in TEM-1 by growing mutant libraries under ampicillin selection — giving us
experimental ground truth for ~5,000 variants, including textbook examples of the tradeoff:
**G238S** extends the substrate range to third-generation cephalosporins but destabilizes the
fold; **M182T** is the canonical compensatory mutation that restores stability without changing
activity.

### Where ESM-2 comes in

[ESM-2][esm2] (Evolutionary Scale Modeling, v2) is a family of protein **language models**
released by Meta AI in 2022–2023. Architecturally, it is a transformer encoder — the same kind
of model as BERT — but trained on protein sequences instead of natural language. The training
objective is **masked language modeling**: random residues are hidden, and the model learns to
reconstruct them from surrounding context. The training corpus is ~65 million unique sequences
from UniRef50, spanning the full diversity of life.

Because the model has to predict masked residues across this corpus, it implicitly learns:

- which residues are **conserved** in a given structural or functional context,
- which **substitutions** are evolutionarily acceptable at a given position,
- and longer-range dependencies between positions (co-evolution, contacts, motifs).

In other words, ESM-2's predicted distribution at a masked position is a learned **evolutionary
prior** — a compressed summary of the substitutions that nature has tolerated in similar
contexts across billions of years. ESM-2 ships in several sizes, from 8M to 15B parameters
(the largest, ESM-2 3B and 15B, were used in ESMFold for structure prediction). For this
tutorial we use the smallest practical model, **ESM-2 35M** (`esm2_t12_35M_UR50D`), which fits on
a CPU laptop; the 650M variant gives sharper predictions if you have a GPU.

The key idea we exploit is this: if a mutation is destabilizing, function-disrupting, or simply
unprecedented in evolution, ESM-2 should assign it a **lower likelihood** than the wild-type
residue. We never train it on TEM-1 fitness data; we just read off the probabilities it learned
during pretraining. This is the **zero-shot** setting.

**Our central question:** Can ESM-2, a protein language model trained only on sequences, predict
which mutations preserve, break, or enhance TEM-1 function — and how does its signal relate to
the stability–activity tradeoff above?

[tokuriki]: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000002
[chen]: https://pmc.ncbi.nlm.nih.gov/articles/PMC1360657/
[studer]: https://www.pnas.org/doi/10.1073/pnas.1310811111
[esm2]: https://www.science.org/doi/10.1126/science.ade2574

## What you'll learn

By the end of this tutorial, you'll have:

- A working pipeline for mutation effect prediction with ESM-2
- An honest sense of where ESM-2 helps and where it doesn't
- A template you can adapt to your own protein of interest

## Prerequisites

- Familiarity with Python and PyTorch
- Basic understanding of transformers (helpful but not required)
- A laptop with 8GB+ RAM, **or** a free Google Colab account

Runtime: ~10 minutes on a CPU laptop.

## Tutorial structure

### Part 1 — Zero-shot scoring with ESM-2

**Notebook:** [`01_zero_shot_scoring.ipynb`](01_zero_shot_scoring.ipynb)

We start without any training. ESM-2 was pretrained on ~65 million protein sequences using
masked language modeling, learning a rich prior over plausible protein sequences. We can exploit
this directly to score mutations — a setup known as
[zero-shot learning](https://en.wikipedia.org/wiki/Zero-shot_learning).

You'll learn to:

- Load ESM-2 and explore its outputs
- Score mutations using **masked marginal likelihood**
- Compute a **per-position constraint profile** that reveals the active site without supervision
- Compare ESM-2 predictions to experimental fitness scores from Firnberg et al.

This part runs comfortably on any laptop (even an 8GB MacBook Air) in ~10 minutes.

**Key takeaway:** ESM-2 already "knows" a surprising amount about protein function from
pretraining alone. But its scores are relative, not calibrated to specific experimental conditions.

## How to run

### Option A: Local installation

Instructions to install dependencies are [here](../../INSTALL.md)

```bash
cd evosite3d/tutorials/esm2-mutation-effects
uv sync
uv run jupyter lab
```

Then open `01_zero_shot_scoring.ipynb`.

### Option B: Google Colab

Click the "Open in Colab" badge below. The notebook installs its dependencies
automatically. A GPU is not required for Part 1 — the 35M model runs on CPU in a few minutes —
but if you want to switch to the 650M model, enable one via
**Runtime → Change runtime type → T4 GPU**.

[![Open Part 1 in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/romainstuder/evosite3d/blob/main/tutorials/esm2-mutation-effects/01_zero_shot_scoring.ipynb) **Part 1: Zero-shot scoring**

## Data

The notebook downloads the Firnberg et al. (2014) deep mutational scanning data from ProteinGym
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

> A note: several of the biology references below are from 2005–2014 — they pre-date ESM-2 by
> a decade. That's because they come out of my own postdoc work on the stability–activity
> tradeoff, and they remain the clearest framing of the problem I know. The biology hasn't
> moved much since; what has moved is our ability to predict its consequences from sequence
> alone, which is the modern half of this tutorial.

- **ESM-2:** Lin et al. (2023). _Evolutionary-scale prediction of atomic-level protein structure._ Science 379(6637):1123-1130.
- **TEM-1 DMS data:** Firnberg et al. (2014). _A comprehensive, high-resolution map of a gene's fitness landscape._ Mol Biol Evol 31(6):1581-1592.
- **Stability–function tradeoff:** Tokuriki, Stricher, Serrano & Tawfik (2008). _How protein stability and new functions trade off._ PLoS Comput Biol 4(2):e1000002.
- **β-lactamase mechanism:** Chen, Shoichet & Bonnet (2005). _Structure, function, and inhibition along the reaction coordinate of CTX-M β-lactamases._ J Am Chem Soc 127(15):5423-5434.
- **RuBisCO motivation:** Studer et al. (2014). _Stability-activity tradeoffs constrain the adaptive evolution of RubisCO._ PNAS 111(6):2223-2228.
- **ProteinGym benchmark:** Notin et al. (2023). _ProteinGym: large-scale benchmarks for protein fitness prediction._ NeurIPS Datasets and Benchmarks.

## License

MIT

## Citation

If this tutorial helps your research or teaching, please cite or link back. Issues and pull
requests welcome.
