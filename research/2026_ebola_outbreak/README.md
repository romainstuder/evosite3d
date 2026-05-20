# Ebola Virus Glycoprotein Interface Mutations: Zaire vs Bundibugyo Strains

## Overview

Comparative analysis of Ebola virus glycoprotein (GP) sequences between Zaire and Bundibugyo strains reveals critical mutations at antibody-binding interface regions that may impact immune recognition and vaccine efficacy.

## Key Findings

### Interface Mutation Landscape

- **21 interface positions** identified using PDB structure 3CSY (GP-antibody complex)
- **8 positions mutated** between strains (38% mutation rate)
- **13 positions conserved** (62% conservation rate)

### Critical Mutations at Interface

The following mutations occur at antibody-binding interface positions:

| Position | Zaire | Bundibugyo | Region            | Impact                |
| -------- | ----- | ---------- | ----------------- | --------------------- |
| 41       | Ser   | Asn        | N-terminal        | Polar→Polar           |
| 503      | Ala   | Ile        | GP1-GP2 junction  | Hydrophobic expansion |
| 504      | Ile   | Thr        | GP1-GP2 junction  | Hydrophobic→Polar     |
| 505      | Val   | Leu        | GP1-GP2 junction  | Hydrophobic shift     |
| 506      | Asn   | Ser        | GP1-GP2 junction  | Polar maintained      |
| 507      | Ala   | Thr        | GP1-GP2 junction  | Small→Polar           |
| 509      | Pro   | Ala        | GP1-GP2 junction  | Conformational change |
| 552      | Asp   | Asn        | C-terminal region | Negative→Neutral      |

![Sequence Alignment](gp_aln.png)
_Figure: Multiple sequence alignment of Ebola virus glycoprotein sequences showing interface mutations between Zaire and Bundibugyo strains. Interface residues are highlighted with mutation status indicated._

### Structural Conservation

**Disulfide Bridge Maintained**: Cysteines at positions **511** and **556** are completely conserved (C511C, C556C), indicating the critical disulfide bridge between these residues remains intact across both strains. This structural element is essential for GP stability and proper folding.

### Hotspot Analysis

The **GP1-GP2 junction region (positions 503-509)** shows the highest mutation density with 6/7 positions altered, suggesting this region may be under immune pressure and could represent an immune escape hotspot.

## Implications

### Vaccine Development

- Current Zaire-based vaccines (rVSV-ZEBOV) may have reduced efficacy against Bundibugyo variants
- The conserved disulfide bridge (C511-C556) represents a stable structural target
- Interface mutations may affect neutralizing antibody binding

### Therapeutic Targeting

- Monoclonal antibodies targeting the GP1-GP2 junction should be evaluated for cross-reactivity
- Conserved interface positions (T42, Q508, K510, C511, P513, N514, H549-L554, C556) represent promising universal targets

### Surveillance Priority

- Monitor GP1-GP2 junction mutations in circulating strains
- Track conservation of structural cysteines as indicators of viral fitness constraints

## Methods

- Interface residues identified using PISA analysis of PDB 3CSY
- Sequence alignment performed with MAFFT
- 5Å distance cutoff for antibody-antigen contacts
- Analysis covers 21 interface positions across GP sequence

## Data Availability

- Alignment files: `gp_aln.fasta`
- Jalview annotation: `interface_annotation.jva`
- Analysis notebook: `identify_mutations.ipynb`

---

_Analysis performed on 2026-05-20 using EvoSite3D computational pipeline_
