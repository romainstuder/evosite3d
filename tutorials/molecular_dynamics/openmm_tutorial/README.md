# OpenMM Ubiquitin Molecular Dynamics Tutorial

A complete Jupyter notebook tutorial for running a 5 ns molecular dynamics (MD) simulation of ubiquitin using OpenMM, plus VMD scripts for visualisation and movie rendering. Runs on CUDA, OpenCL, or CPU — the notebook auto-detects the fastest available platform.

## Table of Contents

- [About Ubiquitin](#about-ubiquitin)
- [What is Molecular Dynamics?](#what-is-molecular-dynamics)
- [Installation](#installation)
- [Workflow Overview](#workflow-overview)
- [Step-by-Step Guide](#step-by-step-guide)
- [Expected Results](#expected-results)
- [Trajectory Visualisation (VMD)](#trajectory-visualisation-vmd)
- [Troubleshooting](#troubleshooting)

---

## About Ubiquitin

**Ubiquitin** is a small (76 amino acids, ~8.5 kDa), highly conserved regulatory protein found in nearly all eukaryotic cells. Its name comes from its "ubiquitous" presence throughout living organisms.

### Why ubiquitin?

- **Small and fast** - Ideal for tutorials and testing (simulates quickly)
- **Well-characterized** - Extensively studied with abundant experimental data
- **Biologically important** - Tags proteins for degradation by the proteasome
- **Stable structure** - Has a compact β-grasp fold that's well-behaved in simulations
- **PDB structure** - High-resolution X-ray structure available (PDB: 1UBQ, 1.8 Å resolution)

### Biological function

Ubiquitin attaches to other proteins (ubiquitination) to mark them for:

- Degradation by the proteasome
- Cellular signaling
- DNA repair
- Cell cycle regulation

This process is so important that the discovery of ubiquitin-mediated protein degradation won the **2004 Nobel Prize in Chemistry**.

---

## What is Molecular Dynamics?

**Molecular Dynamics (MD)** is a computational technique that simulates the physical movements of atoms over time by solving Newton's equations of motion.

### How it works

1. Each atom is treated as a particle with mass, charge, and position
2. Forces between atoms are calculated using a **force field** (mathematical description of bonded and non-bonded interactions)
3. Newton's equations (F = ma) are integrated over tiny time steps (typically 2 femtoseconds)
4. Atomic positions are updated to generate a trajectory through time

### Why do MD?

MD provides **atomic-level insight** into protein behavior that experiments cannot directly observe:

- 🔬 **Protein dynamics** - How proteins move and flex
- 💊 **Drug discovery** - How drugs bind to targets
- 🧬 **Protein folding** - How proteins acquire their 3D structure
- ⚗️ **Conformational changes** - Transitions between states
- 🌡️ **Thermodynamic properties** - Free energies, binding affinities
- 🔗 **Protein-protein interactions** - How proteins interact

### Limitations

- Computationally expensive (typically nanoseconds to microseconds of simulation time)
- Limited accuracy (depends on force field quality)
- Sampling problem (rare events may not be captured in short simulations)

---

## Installation

### Requirements

- Python 3.10+
- OpenMM 8.0+
- MDTraj
- NumPy, Matplotlib
- Jupyter

### Quick Setup

This project uses [uv](https://docs.astral.sh/uv/) for environment and dependency management.

```bash
# From the repo root
uv sync                 # installs everything in pyproject.toml + uv.lock
uv run jupyter lab tutorials/molecular_dynamics/openmm_tutorial/OpenMM_Ubiquitin_Tutorial.ipynb
```

To add the dependencies into a fresh project:

```bash
uv add openmm mdtraj matplotlib numpy jupyter
```

### Optional: Metal GPU on Apple Silicon

OpenMM's Metal platform is a community plugin and is **not** bundled. Install it if you want GPU acceleration on Apple Silicon:

```bash
uv add openmm-metal
```

Otherwise the notebook will fall back to OpenCL (works on most Macs) or CPU.

---

## Workflow Overview

| Step | Process                    | Purpose                                        | Duration (OpenCL on M1) |
| ---- | -------------------------- | ---------------------------------------------- | ----------------------- |
| 1    | Load PDB                   | Get initial structure                          | <1 s                    |
| 2    | System Setup               | Add water, ions, hydrogens                     | ~10 s                   |
| 3    | Energy Minimization        | Fix bad contacts                               | ~15 s                   |
| 4    | NVT Equilibration (100 ps) | Heat to 300 K                                  | ~75 s                   |
| 5    | NPT Equilibration (100 ps) | Equilibrate pressure                           | ~80 s                   |
| 6    | **Production MD (5 ns)**   | **Sample equilibrium dynamics**                | **~60 min**             |
| 7    | Analysis                   | Concatenate trajectories, compute RMSD/RMSF/Rg | ~10 s                   |
| 8    | Visualization              | Generate stage-annotated plots                 | <5 s                    |

**Total time:** ~65 min on Apple OpenCL. CUDA is comparable; pure CPU is 5-10× slower. Set `sim_ns = 1` in cell 21 to fall back to a 1 ns / ~15 min run for quick iteration.

---

## Step-by-Step Guide

### Step 1: Load Structure

**What:** Load the PDB file 1UBQ from the Protein Data Bank

**Why:** This is our starting point - the experimentally-determined 3D structure of ubiquitin. The crystal structure provides initial atomic coordinates for all 76 residues.

### Step 2: System Setup

**What:**

- Add water molecules around the protein (TIP3P water model)
- Add ions (Na+ and Cl- at physiological 0.15 M concentration)
- Add hydrogens (PDB structures often lack hydrogens)
- Define the simulation box with periodic boundary conditions

**Why:**

- 💧 **Water** - Proteins exist in aqueous environment; water dramatically affects their behavior
- 🧂 **Ions** - Physiological salt concentration mimics cellular conditions
- ⚛️ **Hydrogens** - Essential for accurate bonding, especially hydrogen bonds
- 📦 **Periodic box** - Simulates an infinite system without edge effects

### Step 3: Energy Minimization

**What:** Use steepest descent algorithm to find a low-energy structure

**Why:** The initial structure often has small atomic clashes (steric overlaps) from:

- Added hydrogens being in default positions
- Water molecules placed in unfavorable orientations
- Crystal packing artifacts from the experimental structure

Minimization removes these clashes by moving atoms downhill on the energy surface. Without this step, the simulation would crash or behave erratically.

### Step 4: NVT Equilibration (100 ps)

**What:** Gradually heat the system to 300 K (room temperature) at **constant volume**

**Why:**

- 🌡️ **Temperature ramp** - Sudden heating causes structural damage
- ⚖️ **Constant volume (NVT)** - Easier to control initial heating
- 🎯 **Goal** - Reach target temperature while maintaining structure

This step uses a **Langevin thermostat** to control temperature, which adds friction and random forces to maintain the target temperature.

### Step 5: NPT Equilibration (100 ps)

**What:** Equilibrate at **constant pressure** (1 bar) and temperature (300 K)

**Why:**

- 🎈 **Pressure equilibration** - Allows box size to adjust to correct density (~1000 kg/m³)
- 🏠 **Physiological conditions** - Atmospheric pressure matches cellular environment
- 📊 **Density convergence** - Ensures proper water density before production

A **Monte Carlo barostat** randomly adjusts the box size to maintain target pressure.

### Step 6: Production MD (5 ns)

**What:** Run the actual simulation, collecting trajectory data

**Why:** This is where the science happens. Once equilibrated, the system explores conformational space, and we record:

- 🎬 **Trajectory** (DCD file) - Atomic positions over time, sampled every 10 ps
- ⚡ **Energy** (text file) - Potential, kinetic, total energy
- 🌡️ **Temperature & density** - System properties
- 💾 **Checkpoint** - Periodic snapshot for crash recovery (`CheckpointReporter`)

The simulation length is controlled by the `sim_ns` variable in cell 21. 5 ns gives a clear equilibrium plateau on the RMSD plot; 1 ns is enough for a first run. For large conformational changes you'd need μs-ms simulations, beyond the scope of this tutorial.

### Step 7: Trajectory Analysis

NVT, NPT, and Production trajectories are joined into a single timeline before analysis (`md.join`). Periodic-boundary jumps are corrected with `traj.image_molecules`, and all metrics are computed on the **protein only** (with backbone superposition) — never on the bulk water, which would render RMSD/RMSF meaningless.

We compute three key metrics:

#### RMSD (Root Mean Square Deviation)

**What it measures:** How much the backbone deviates from the starting point (after superposition)

**Interpretation:**

- Low RMSD (~0.10–0.15 nm for ubiquitin) = Structure is stable
- Rising RMSD = Structure is changing
- Plateau = Equilibrium reached

#### RMSF (Root Mean Square Fluctuation)

**What it measures:** Per-residue flexibility (how much each amino acid moves)

**Interpretation:**

- Low RMSF = Rigid regions (typically core, β-sheets)
- High RMSF = Flexible regions (typically loops, termini)
- Tells us which parts of the protein are dynamic

#### Radius of Gyration (Rg)

**What it measures:** Overall protein compactness

**Interpretation:**

- Stable Rg = Protein maintains its fold
- Increasing Rg = Protein is unfolding/expanding
- Decreasing Rg = Protein is compacting

### Step 8: Visualization

Generate publication-quality plots showing:

- RMSD vs time (stability) — with stage shading (NVT / NPT / Production)
- RMSF per residue (flexibility map)
- Rg vs time (compactness) — with stage shading
- Summary statistics

---

## Expected Results

After running the complete notebook with ubiquitin:

| Metric              | Expected Value         | Meaning                          |
| ------------------- | ---------------------- | -------------------------------- |
| Final backbone RMSD | 0.10-0.15 nm           | Protein is stable                |
| Mean Rg             | ~1.18 nm               | Compact β-grasp fold             |
| Max RMSF            | 0.3-0.4 nm @ res 72-76 | Flexible C-terminal Gly-Gly tail |
| Temperature         | 300 ± 5 K              | Well-equilibrated                |
| Density             | ~1.01 g/mL             | Realistic water density          |
| Performance         | ~120 ns/day            | Apple OpenCL (M1)                |

### Biological Insights

You should observe:

- 🔒 **Stable core** - β-sheets show low RMSF
- 🌊 **Flexible C-terminus** - The Gly-Gly tail (where ubiquitin attaches to substrates) is highly mobile
- 🔄 **Loop dynamics** - Loops between secondary structure elements fluctuate
- 💪 **Compact structure** - Ubiquitin maintains its characteristic fold

---

## Trajectory Visualisation (VMD)

Two helper Tcl scripts wrap a sensible default view and a turntable movie pipeline.

### `setup_view.tcl` — interactive scene

```bash
vmd -e setup_view.tcl
```

What it does:

- Loads `output/step3_minimised.pdb` + `output/step6_trajectory.dcd`
- Hides solvent (only `protein` selections get reps)
- Builds a NewCartoon for the backbone
- Aligns every frame to frame 0 on the backbone, so the protein doesn't drift
- **Colours the cartoon by per-residue RMSF** using `output/step7_rmsf.txt` (BGYR scale: rigid β-sheet core in blue, flexible C-terminal tail in red)
- White background, orthographic projection, axes off

Configurable variables at the top of the script let you repoint the files or set `RMSF_FILE ""` to fall back to standard secondary-structure coloring.

### `make_movie.tcl` — render a turntable movie

```bash
vmd -e make_movie.tcl
```

Sources `setup_view.tcl` first (so you get the same RMSF-coloured scene), then renders one TGA per frame while doing a 360° turntable rotation. Configurable knobs at the top:

| Variable           | Default | Purpose                                       |
| ------------------ | ------- | --------------------------------------------- |
| `STRIDE`           | 1       | Render every Nth trajectory frame             |
| `ROTATE_TOTAL_DEG` | 360.0   | Total camera rotation across the movie        |
| `ROTATE_AXIS`      | `"y"`   | `x` (tumble), `y` (turntable), `z` (roll)     |
| `ZOOM`             | 2.0     | View scale; >1 zooms in to crop empty solvent |

Frames land in `output/movie_frames/`. Stitch into MP4 with ffmpeg:

```bash
ffmpeg -framerate 25 -i output/movie_frames/frame_%04d.tga \
       -c:v libx264 -pix_fmt yuv420p output/movie.mp4
```

The script prints this command at the end as a reminder.

---

## Force Field Details

**AMBER99SB-ILDN** is used because:

- ✅ Industry-standard for protein simulations
- ✅ Well-validated against NMR data
- ✅ Accurate for backbone dynamics
- ✅ Improved side chain rotamers (ILDN modifications)

**TIP3P water** is used because:

- ✅ Computationally efficient
- ✅ Compatible with AMBER force fields
- ✅ Standard choice for biomolecular MD

---

## Troubleshooting

### `amber99sbildn.xml not found`

Force field XMLs ship with OpenMM itself — if missing, your install is broken. Reinstall:

```bash
uv sync --reinstall openmm
```

### Platform selection

After running the platform-detection cell, the notebook prints which platform it chose. Order of preference is CUDA → Metal → OpenCL → CPU. To check at any later point:

```python
print(simulation.context.getPlatform().getName())
```

### `OpenMMException: No compatible OpenCL platform is available`

Misleading error — OpenCL is usually fine, but Apple's driver does **not** accept `{'Precision': 'mixed'}`. The notebook only sets that property on CUDA. If you've customised it, restrict mixed precision to CUDA only.

### Simulation crashes (NaN forces, infinite energy)

- Ensure energy minimization completed successfully (`Final energy` should be a large negative number, not NaN)
- Check that hydrogens were added before `createSystem`
- Verify all residues have force field parameters (non-standard residues will fail)
- Reduce timestep from 2 fs to 1 fs while diagnosing

### Out of memory

Reduce solvent padding (`padding=0.8*nanometer` instead of `1.0`) or force single precision on CUDA:

```python
properties = {'Precision': 'single'}
```

---

## Next Steps

After completing this tutorial:

1. **Run longer** - Bump `sim_ns` from 5 to 50 ns (~10 h on OpenCL) for converged sampling of loop dynamics
2. **Run multiple replicas** - Statistical significance requires 3-5 independent simulations with different random seeds
3. **Different proteins** - Apply this workflow to your protein of interest
4. **Advanced analysis** - DSSP (secondary structure), hydrogen bonds, contact maps
5. **Enhanced sampling** - Replica exchange MD, metadynamics, umbrella sampling
6. **Free energy calculations** - Binding affinities, conformational free energies

---

## References

### Software

- **OpenMM**: Eastman et al. (2017) _PLOS Comp Biol_ 13:e1005659
- **MDTraj**: McGibbon et al. (2015) _Biophys J_ 109:1528-1532
- **AMBER99SB-ILDN**: Lindorff-Larsen et al. (2010) _Proteins_ 78:1950-1958

### Ubiquitin

- **Structure**: Vijay-Kumar et al. (1987) _J Mol Biol_ 194:531-544 (PDB: 1UBQ)
- **Nobel Prize 2004**: Aaron Ciechanover, Avram Hershko, Irwin Rose

### Tutorials

- OpenMM User Guide: http://docs.openmm.org
- MDTraj Documentation: http://mdtraj.org

---

## License

This tutorial is provided for educational purposes. Feel free to adapt and share.

## Contact

For questions or issues, open an issue on GitHub or contact the author.

---

**Happy simulating!** 🧬⚛️
