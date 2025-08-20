# Comprehensive Guide to Protein Structure Quality Assessment

## Table of Contents

1. [Introduction](#introduction)
2. [Resolution-Based Metrics](#resolution-based-metrics)
3. [Geometric Validation](#geometric-validation)
4. [Model vs. Data Agreement](#model-vs-data-agreement)
5. [Structural Completeness](#structural-completeness)
6. [B-factor Analysis](#b-factor-analysis)
7. [Comparative Assessment](#comparative-assessment)
8. [Practical Implementation](#practical-implementation)
9. [Automated Quality Assessment](#automated-quality-assessment)
10. [Best Practices and Interpretation](#best-practices-and-interpretation)

---

## Introduction

Protein structure quality assessment is crucial for determining the reliability of structural data
for downstream applications. Quality metrics fall into several categories: experimental data
quality, geometric validity, completeness, and agreement with known structural principles.

**File Formats:**

- **PDB format**: Legacy text format, limited precision
- **mmCIF format**: Modern format with higher precision and extensibility
- **Both contain**: Atomic coordinates, experimental metadata, validation metrics

---

## Resolution-Based Metrics

### 1. Resolution

**Description**: Measure of the finest detail visible in the electron density map (X-ray) or
reconstruction (cryo-EM).

**Interpretation**:

- **X-ray crystallography**: < 2.0 Å (excellent), 2.0-3.0 Å (good), > 3.0 Å (lower detail)
- **Cryo-EM**: < 3.0 Å (excellent), 3.0-4.0 Å (good), > 4.0 Å (lower detail)
- **NMR**: Not applicable (reports ensemble of structures)

**Assessment Tools**:

```bash
# Extract resolution from PDB header
grep "^REMARK   2 RESOLUTION" structure.pdb

# Using BioPython
from Bio.PDB import PDBParser
parser = PDBParser()
structure = parser.get_structure("protein", "structure.pdb")
resolution = structure.header.get('resolution')
print(f"Resolution: {resolution} Å")
```

### 2. Wilson B-factor

**Description**: Overall temperature factor indicating atomic displacement parameters and data
quality.

**Interpretation**:

- **Low values (< 20 Å²)**: Well-ordered, high-quality crystals
- **High values (> 50 Å²)**: Disorder, radiation damage, or poor crystals

```python
# Extract Wilson B-factor from mmCIF
import gemmi
structure = gemmi.read_structure("structure.cif")
wilson_b = structure.info.get("_reflns.B_iso_Wilson_estimate")
```

---

## Geometric Validation

### 3. Ramachandran Plot Analysis

**Description**: φ/ψ dihedral angle distribution validation against known protein geometry.

**Metrics**:

- **Favored regions**: > 95% expected
- **Allowed regions**: 4-5% acceptable
- **Outliers**: < 1% for high-quality structures

**Implementation**:

```python
# Using BioPython
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.Polypeptide import protein_letters_3to1
import numpy as np

def calculate_dihedral(p1, p2, p3, p4):
    """Calculate dihedral angle between four points"""
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    n1 = n1 / np.linalg.norm(n1)
    n2 = n2 / np.linalg.norm(n2)

    m1 = np.cross(n1, b2/np.linalg.norm(b2))

    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    return np.degrees(np.arctan2(y, x))

def ramachandran_analysis(structure_file):
    parser = PDBParser()
    structure = parser.get_structure("protein", structure_file)

    phi_psi_angles = []

    for model in structure:
        for chain in model:
            residues = list(chain.get_residues())
            for i in range(1, len(residues)-1):
                try:
                    # Get atoms for phi and psi calculation
                    prev_c = residues[i-1]['C'].coord
                    curr_n = residues[i]['N'].coord
                    curr_ca = residues[i]['CA'].coord
                    curr_c = residues[i]['C'].coord
                    next_n = residues[i+1]['N'].coord

                    # Calculate phi and psi angles
                    phi = calculate_dihedral(prev_c, curr_n, curr_ca, curr_c)
                    psi = calculate_dihedral(curr_n, curr_ca, curr_c, next_n)

                    phi_psi_angles.append((phi, psi, residues[i].resname))
                except KeyError:
                    continue

    return phi_psi_angles

# Usage
angles = ramachandran_analysis("structure.pdb")
```

**Using MolProbity via command line**:

```bash
# Install MolProbity suite
phenix.molprobity structure.pdb

# Or use online MolProbity server
wget http://molprobity.biochem.duke.edu/
```

### 4. Bond Lengths and Angles

**Description**: Validation of covalent geometry against ideal values.

**Expected Values**:

- **Bond length RMSD**: < 0.02 Å for high-quality structures
- **Bond angle RMSD**: < 2.0° for high-quality structures

```python
# Using Gemmi for geometric validation
import gemmi

def validate_geometry(structure_file):
    structure = gemmi.read_structure(structure_file)

    # Get restraints (ideal geometry)
    try:
        restraints = gemmi.read_monomer_lib()

        bond_deviations = []
        angle_deviations = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    # Check bond lengths
                    for atom1, atom2 in get_bonds(residue):
                        if atom1 and atom2:
                            distance = atom1.pos.dist(atom2.pos)
                            ideal_dist = get_ideal_bond_length(atom1.name, atom2.name)
                            if ideal_dist:
                                deviation = abs(distance - ideal_dist)
                                bond_deviations.append(deviation)

        return {
            'bond_rmsd': np.sqrt(np.mean([d**2 for d in bond_deviations])),
            'angle_rmsd': np.sqrt(np.mean([d**2 for d in angle_deviations]))
        }
    except:
        print("Restraints not available")
        return None
```

### 5. Clash Detection

**Description**: Identification of atoms that are too close together (van der Waals overlaps).

**Criteria**:

- **Serious clashes**: > 0.4 Å overlap
- **Acceptable**: < 1% of atoms involved in serious clashes

```python
def detect_clashes(structure_file, cutoff=0.4):
    parser = PDBParser()
    structure = parser.get_structure("protein", structure_file)

    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element != 'H':  # Skip hydrogens
                        atoms.append(atom)

    clashes = []
    vdw_radii = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80}  # Simplified

    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            atom1, atom2 = atoms[i], atoms[j]

            # Skip bonded atoms (simplified check)
            if (atom1.parent == atom2.parent and
                abs(atom1.serial_number - atom2.serial_number) <= 2):
                continue

            distance = atom1 - atom2
            vdw_sum = vdw_radii.get(atom1.element, 1.7) + vdw_radii.get(atom2.element, 1.7)

            overlap = vdw_sum - distance
            if overlap > cutoff:
                clashes.append({
                    'atom1': atom1,
                    'atom2': atom2,
                    'overlap': overlap,
                    'distance': distance
                })

    return clashes
```

---

## Model vs. Data Agreement

### 6. R-factors (X-ray crystallography)

**Description**: Measure of agreement between observed and calculated structure factors.

**Types**:

- **R-work**: Calculated from working set (used in refinement)
- **R-free**: Calculated from test set (not used in refinement)
- **R-free - R-work**: Should be < 5% for valid refinement

**Interpretation**:

- **Excellent**: R-free < 20%, gap < 3%
- **Good**: R-free 20-25%, gap < 5%
- **Questionable**: R-free > 30% or gap > 7%

```bash
# Extract R-factors from PDB
grep "^REMARK   3   R VALUE" structure.pdb
grep "^REMARK   3   FREE R VALUE" structure.pdb

# Using BioPython
def extract_r_factors(pdb_file):
    with open(pdb_file, 'r') as f:
        for line in f:
            if "R VALUE     (WORKING + TEST SET)" in line:
                r_work = float(line.split()[-1])
            elif "FREE R VALUE" in line:
                r_free = float(line.split()[-1])

    return {'r_work': r_work, 'r_free': r_free, 'gap': r_free - r_work}
```

### 7. Cross-Correlation (Cryo-EM)

**Description**: Measure of agreement between atomic model and electron density map.

**Metrics**:

- **FSC (Fourier Shell Correlation)**: Overall model-map correlation
- **Local correlation**: Per-residue agreement
- **Q-score**: Quality measure at atomic level

```python
# Using ChimeraX command line for cryo-EM validation
# chimeraX --nogui --cmd "open structure.pdb; open map.mrc; measure correlation #1 #2"

# Or using PHENIX
# phenix.real_space_refine model.pdb map.ccp4
```

---

## Structural Completeness

### 8. Missing Residues and Atoms

**Description**: Assessment of model completeness relative to the full sequence.

```python
def assess_completeness(structure_file, sequence=None):
    parser = PDBParser()
    structure = parser.get_structure("protein", structure_file)

    results = {}

    for model in structure:
        for chain in model:
            residues = list(chain.get_residues())

            # Count missing atoms per residue
            missing_atoms = 0
            total_expected = 0

            for residue in residues:
                expected_atoms = get_expected_atoms(residue.resname)
                present_atoms = len(residue)

                missing_atoms += max(0, expected_atoms - present_atoms)
                total_expected += expected_atoms

            # Check for missing residues (gaps in numbering)
            residue_numbers = [res.id[1] for res in residues]
            missing_residues = []

            if len(residue_numbers) > 1:
                for i in range(min(residue_numbers), max(residue_numbers) + 1):
                    if i not in residue_numbers:
                        missing_residues.append(i)

            results[chain.id] = {
                'missing_atoms': missing_atoms,
                'completeness': (total_expected - missing_atoms) / total_expected * 100,
                'missing_residues': missing_residues
            }

    return results

def get_expected_atoms(resname):
    # Simplified - should use actual residue definitions
    standard_residues = {
        'ALA': 5, 'CYS': 6, 'ASP': 8, 'GLU': 9, 'PHE': 11,
        'GLY': 4, 'HIS': 10, 'ILE': 8, 'LYS': 9, 'LEU': 8,
        'MET': 8, 'ASN': 8, 'PRO': 7, 'GLN': 9, 'ARG': 11,
        'SER': 6, 'THR': 7, 'VAL': 7, 'TRP': 14, 'TYR': 12
    }
    return standard_residues.get(resname, 5)  # Default estimate
```

### 9. Occupancy Analysis

**Description**: Assessment of atom occupancy values and alternate conformations.

```python
def analyze_occupancy(structure_file):
    parser = PDBParser()
    structure = parser.get_structure("protein", structure_file)

    occupancy_stats = {
        'low_occupancy': [],  # < 0.5
        'zero_occupancy': [],  # = 0.0
        'alternate_conformations': []
    }

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.occupancy < 0.5:
                        occupancy_stats['low_occupancy'].append(atom)
                    if atom.occupancy == 0.0:
                        occupancy_stats['zero_occupancy'].append(atom)
                    if atom.altloc not in ['', 'A']:
                        occupancy_stats['alternate_conformations'].append(atom)

    return occupancy_stats
```

---

## B-factor Analysis

### 10. Temperature Factor Distribution

**Description**: Analysis of atomic displacement parameters indicating flexibility and disorder.

```python
import matplotlib.pyplot as plt
import numpy as np

def analyze_b_factors(structure_file):
    parser = PDBParser()
    structure = parser.get_structure("protein", structure_file)

    b_factors = []
    ca_b_factors = []  # C-alpha atoms only

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    b_factors.append(atom.bfactor)
                    if atom.name == 'CA':
                        ca_b_factors.append((residue.id[1], atom.bfactor))

    stats = {
        'mean_b': np.mean(b_factors),
        'median_b': np.median(b_factors),
        'std_b': np.std(b_factors),
        'min_b': np.min(b_factors),
        'max_b': np.max(b_factors),
        'ca_profile': ca_b_factors
    }

    # Plot B-factor profile
    plt.figure(figsize=(12, 6))

    # Histogram
    plt.subplot(1, 2, 1)
    plt.hist(b_factors, bins=50, alpha=0.7)
    plt.xlabel('B-factor (Ų)')
    plt.ylabel('Frequency')
    plt.title('B-factor Distribution')

    # C-alpha profile
    plt.subplot(1, 2, 2)
    if ca_b_factors:
        residues, b_vals = zip(*sorted(ca_b_factors))
        plt.plot(residues, b_vals, 'b-', linewidth=1)
        plt.xlabel('Residue Number')
        plt.ylabel('C-alpha B-factor (Ų)')
        plt.title('B-factor Profile Along Chain')

    plt.tight_layout()
    plt.savefig('b_factor_analysis.png', dpi=300)

    return stats
```

---

## Comparative Assessment

### 11. Structure Similarity Metrics

**Description**: Comparison with homologous structures for validation.

```python
# Using PyMOL for structure alignment and RMSD calculation
def compare_structures(ref_structure, test_structure):
    """
    PyMOL script for structure comparison
    """
    pymol_script = f"""
    load {ref_structure}, ref
    load {test_structure}, test

    # Align structures
    align test, ref

    # Calculate RMSD
    rms_cur test, ref

    # Save alignment statistics
    save alignment_results.pse
    """

    with open('compare_structures.pml', 'w') as f:
        f.write(pymol_script)

    # Run PyMOL script
    import subprocess
    subprocess.run(['pymol', '-c', 'compare_structures.pml'])

# Using BioPython for RMSD calculation
from Bio.PDB import Superimposer

def calculate_rmsd(structure1_file, structure2_file):
    parser = PDBParser()

    struct1 = parser.get_structure("ref", structure1_file)
    struct2 = parser.get_structure("test", structure2_file)

    # Extract C-alpha atoms
    atoms1 = []
    atoms2 = []

    residues1 = list(struct1[0]['A'].get_residues())
    residues2 = list(struct2[0]['A'].get_residues())

    # Match residues by number
    for res1 in residues1:
        for res2 in residues2:
            if res1.id[1] == res2.id[1]:
                if 'CA' in res1 and 'CA' in res2:
                    atoms1.append(res1['CA'])
                    atoms2.append(res2['CA'])
                break

    if len(atoms1) != len(atoms2):
        print("Warning: Different number of matching atoms")

    # Superimpose and calculate RMSD
    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    sup.apply(atoms2)

    return {
        'rmsd': sup.rms,
        'n_atoms': len(atoms1),
        'transformation': sup.rotran
    }
```

---

## Practical Implementation

### 12. Automated Quality Assessment Pipeline

```python
#!/usr/bin/env python3
"""
Comprehensive protein structure quality assessment pipeline
"""

import sys
import os
from pathlib import Path
import json
from Bio.PDB import PDBParser
import numpy as np

class StructureQualityAssessment:
    def __init__(self, structure_file):
        self.structure_file = structure_file
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure("protein", structure_file)
        self.results = {}

    def run_full_assessment(self):
        """Run complete quality assessment pipeline"""
        print("Running comprehensive structure quality assessment...")

        # Basic metadata
        self.results['metadata'] = self.extract_metadata()

        # Geometric validation
        self.results['ramachandran'] = self.ramachandran_analysis()
        self.results['clashes'] = self.detect_clashes()
        self.results['geometry'] = self.validate_geometry()

        # Completeness analysis
        self.results['completeness'] = self.assess_completeness()
        self.results['occupancy'] = self.analyze_occupancy()

        # B-factor analysis
        self.results['b_factors'] = self.analyze_b_factors()

        # Overall quality score
        self.results['quality_score'] = self.calculate_overall_score()

        return self.results

    def extract_metadata(self):
        """Extract basic structure metadata"""
        header = self.structure.header

        metadata = {
            'resolution': header.get('resolution'),
            'method': header.get('structure_method'),
            'r_work': None,
            'r_free': None,
            'deposition_date': header.get('deposition_date'),
            'release_date': header.get('release_date')
        }

        # Extract R-factors from PDB file
        try:
            with open(self.structure_file, 'r') as f:
                for line in f:
                    if "R VALUE     (WORKING + TEST SET)" in line:
                        metadata['r_work'] = float(line.split()[-1])
                    elif "FREE R VALUE" in line:
                        metadata['r_free'] = float(line.split()[-1])
        except:
            pass

        return metadata

    def calculate_overall_score(self):
        """Calculate overall quality score (0-100)"""
        score = 100

        # Resolution penalty
        if self.results['metadata']['resolution']:
            res = self.results['metadata']['resolution']
            if res > 3.0:
                score -= 20
            elif res > 2.5:
                score -= 10

        # R-factor penalty
        if self.results['metadata']['r_free']:
            r_free = self.results['metadata']['r_free']
            if r_free > 0.30:
                score -= 25
            elif r_free > 0.25:
                score -= 15

        # Ramachandran penalty
        rama = self.results['ramachandran']
        if rama['outlier_percent'] > 2:
            score -= 15
        elif rama['outlier_percent'] > 1:
            score -= 5

        # Completeness penalty
        comp = self.results['completeness']
        avg_completeness = np.mean([chain['completeness'] for chain in comp.values()])
        if avg_completeness < 90:
            score -= 20
        elif avg_completeness < 95:
            score -= 10

        return max(0, score)

    def generate_report(self, output_file=None):
        """Generate comprehensive quality report"""
        if not self.results:
            self.run_full_assessment()

        report = []
        report.append("=" * 60)
        report.append("PROTEIN STRUCTURE QUALITY ASSESSMENT REPORT")
        report.append("=" * 60)
        report.append(f"Structure file: {self.structure_file}")
        report.append("")

        # Overall score
        score = self.results['quality_score']
        if score >= 80:
            quality = "EXCELLENT"
        elif score >= 60:
            quality = "GOOD"
        elif score >= 40:
            quality = "FAIR"
        else:
            quality = "POOR"

        report.append(f"OVERALL QUALITY SCORE: {score:.1f}/100 ({quality})")
        report.append("")

        # Metadata section
        report.append("EXPERIMENTAL DATA:")
        report.append("-" * 20)
        meta = self.results['metadata']
        if meta['resolution']:
            report.append(f"Resolution: {meta['resolution']:.2f} Å")
        report.append(f"Method: {meta['method']}")
        if meta['r_work'] and meta['r_free']:
            report.append(f"R-work: {meta['r_work']:.3f}")
            report.append(f"R-free: {meta['r_free']:.3f}")
            report.append(f"R-gap: {meta['r_free'] - meta['r_work']:.3f}")
        report.append("")

        # Ramachandran analysis
        report.append("RAMACHANDRAN ANALYSIS:")
        report.append("-" * 22)
        rama = self.results['ramachandran']
        report.append(f"Favored: {rama['favored_percent']:.1f}%")
        report.append(f"Allowed: {rama['allowed_percent']:.1f}%")
        report.append(f"Outliers: {rama['outlier_percent']:.1f}%")
        report.append("")

        # Completeness
        report.append("MODEL COMPLETENESS:")
        report.append("-" * 19)
        for chain_id, comp in self.results['completeness'].items():
            report.append(f"Chain {chain_id}: {comp['completeness']:.1f}% complete")
            if comp['missing_residues']:
                report.append(f"  Missing residues: {comp['missing_residues'][:5]}{'...' if len(comp['missing_residues']) > 5 else ''}")
        report.append("")

        # B-factor analysis
        report.append("B-FACTOR ANALYSIS:")
        report.append("-" * 18)
        bf = self.results['b_factors']
        report.append(f"Mean B-factor: {bf['mean_b']:.1f} Ų")
        report.append(f"B-factor range: {bf['min_b']:.1f} - {bf['max_b']:.1f} Ų")
        report.append("")

        # Recommendations
        report.append("RECOMMENDATIONS:")
        report.append("-" * 15)
        if score >= 80:
            report.append("✓ Structure is of excellent quality for most applications")
        elif score >= 60:
            report.append("✓ Structure is suitable for most analyses")
            report.append("⚠ Consider resolution limitations for detailed studies")
        else:
            report.append("⚠ Structure has quality issues")
            report.append("⚠ Use with caution for detailed structural analysis")
            report.append("⚠ Consider alternative structures if available")

        report_text = "\n".join(report)

        if output_file:
            with open(output_file, 'w') as f:
                f.write(report_text)
        else:
            print(report_text)

        return report_text

# Command-line usage
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python structure_qa.py <structure_file>")
        sys.exit(1)

    structure_file = sys.argv[1]
    if not os.path.exists(structure_file):
        print(f"Error: File {structure_file} not found")
        sys.exit(1)

    # Run assessment
    qa = StructureQualityAssessment(structure_file)
    results = qa.run_full_assessment()

    # Generate report
    report_file = Path(structure_file).stem + "_quality_report.txt"
    qa.generate_report(report_file)

    # Save detailed results as JSON
    json_file = Path(structure_file).stem + "_quality_results.json"
    with open(json_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    print(f"\nDetailed report saved to: {report_file}")
    print(f"JSON results saved to: {json_file}")
```

### 13. Unix Command Line Tools

```bash
#!/bin/bash
# Structure quality assessment using command-line tools

STRUCTURE_FILE=$1

if [ -z "$STRUCTURE_FILE" ]; then
    echo "Usage: $0 <structure_file>"
    exit 1
fi

echo "=== PROTEIN STRUCTURE QUALITY ASSESSMENT ==="
echo "File: $STRUCTURE_FILE"
echo

# Basic information extraction
echo "=== BASIC INFORMATION ==="
echo "Resolution:"
grep "^REMARK   2 RESOLUTION" "$STRUCTURE_FILE"

echo "R-factors:"
grep "^REMARK   3.*R VALUE" "$STRUCTURE_FILE"

echo "Experimental method:"
grep "^EXPDTA" "$STRUCTURE_FILE"
echo

# Count atoms and residues
echo "=== MODEL STATISTICS ==="
echo "Total atoms: $(grep "^ATOM" "$STRUCTURE_FILE" | wc -l)"
echo "Protein atoms: $(grep "^ATOM" "$STRUCTURE_FILE" | awk '$4!="HOH"' | wc -l)"
echo "Water molecules: $(grep "^ATOM" "$STRUCTURE_FILE" | awk '$4=="HOH"' | wc -l)"
echo "Unique residues: $(grep "^ATOM" "$STRUCTURE_FILE" | awk '{print $4}' | sort -u | wc -l)"
echo

# Check for missing atoms (simplified)
echo "=== COMPLETENESS CHECK ==="
echo "Residues missing backbone atoms:"
grep "^ATOM" "$STRUCTURE_FILE" | awk '{print $4, $6, $3}' | sort | uniq | \
while read res_name res_num atom_name; do
    if [ "$atom_name" = "CA" ]; then
        # Check if N and C atoms exist for this residue
        n_count=$(grep "^ATOM" "$STRUCTURE_FILE" | awk -v rn="$res_num" -v an="N" '$6==rn && $3==an' | wc -l)
        c_count=$(grep "^ATOM" "$STRUCTURE_FILE" | awk -v rn="$res_num" -v an="C" '$6==rn && $3==an' | wc -l)

        if [ $n_count -eq 0 ] || [ $c_count -eq 0 ]; then
            echo "  Residue $res_name $res_num: missing backbone atoms"
        fi
    fi
done
echo

# B-factor analysis
echo "=== B-FACTOR ANALYSIS ==="
b_factors=$(grep "^ATOM" "$STRUCTURE_FILE" | awk '{print $11}')
echo "B-factor statistics:"
echo "$b_factors" | awk '
    BEGIN { min=999; max=0; sum=0; count=0 }
    {
        if ($1 < min) min = $1;
        if ($1 > max) max = $1;
        sum += $1;
        count++;
    }
    END {
        if (count > 0) {
            mean = sum/count;
            print "  Mean: " mean " Ų";
            print "  Range: " min " - " max " Ų";
        }
    }
'
echo

# Check for alternate conformations
echo "=== ALTERNATE CONFORMATIONS ==="
alt_conf=$(grep "^ATOM" "$STRUCTURE_FILE" | awk '$5!=" " && $5!="A"' | wc -l)
if [ $alt_conf -gt 0 ]; then
    echo "Atoms with alternate conformations: $alt_conf"
    echo "Residues with alternates:"
    grep "^ATOM" "$STRUCTURE_FILE" | awk '$5!=" " && $5!="A" {print $4, $6}' | sort -u
else
    echo "No alternate conformations found"
fi
echo

# Check for zero occupancy atoms
echo "=== OCCUPANCY ANALYSIS ==="
zero_occ=$(grep "^ATOM" "$STRUCTURE_FILE" | awk '$10==0.00' | wc -l)
low_occ=$(grep "^ATOM" "$STRUCTURE_FILE" | awk '$10<0.50 && $10>0.00' | wc -l)

echo "Zero occupancy atoms: $zero_occ"
echo "Low occupancy atoms (< 0.5): $low_occ"

if [ $zero_occ -gt 0 ]; then
    echo "Residues with zero occupancy atoms:"
    grep "^ATOM" "$STRUCTURE_FILE" | awk '$10==0.00 {print $4, $6}' | sort -u | head -5
fi
echo

echo "=== QUALITY ASSESSMENT COMPLETE ==="
```

### 14. Advanced Validation with External Tools

```bash
#!/bin/bash
# Advanced structure validation using external tools

STRUCTURE_FILE=$1
BASE_NAME=$(basename "$STRUCTURE_FILE" .pdb)

# Create output directory
mkdir -p "${BASE_NAME}_validation"
cd "${BASE_NAME}_validation"

echo "Running advanced structure validation..."

# 1. MolProbity validation (if available)
if command -v phenix.molprobity &> /dev/null; then
    echo "Running MolProbity validation..."
    phenix.molprobity "../$STRUCTURE_FILE" > molprobity_results.txt 2>&1
    echo "✓ MolProbity results saved to molprobity_results.txt"
fi

# 2. SAVES validation (if available)
if command -v saves &> /dev/null; then
    echo "Running SAVES validation..."
    saves "../$STRUCTURE_FILE" > saves_results.txt 2>&1
    echo "✓ SAVES results saved to saves_results.txt"
fi

# 3. Ramachandran plot with Rampage (if available)
if command -v rampage &> /dev/null; then
    echo "Generating Ramachandran plot..."
    rampage -file "../$STRUCTURE_FILE" -ps ramachandran.ps
    echo "✓ Ramachandran plot saved to ramachandran.ps"
fi

# 4. Structure factor validation (if MTZ file available)
if [ -f "../${BASE_NAME}.mtz" ]; then
    if command -v phenix.model_vs_data &> /dev/null; then
        echo "Running model vs data validation..."
        phenix.model_vs_data "../$STRUCTURE_FILE" "../${BASE_NAME}.mtz" > model_vs_data.txt 2>&1
        echo "✓ Model vs data results saved to model_vs_data.txt"
    fi
fi

# 5. Generate summary report
cat > validation_summary.txt << EOF
STRUCTURE VALIDATION SUMMARY
============================
Structure: $STRUCTURE_FILE
Date: $(date)

Files generated:
$(ls -la *.txt *.ps 2>/dev/null | grep -v "^total")

Recommendations:
1. Check molprobity_results.txt for geometric validation
2. Review Ramachandran plot for backbone conformations
3. Examine B-factor distribution for structural flexibility
4. Validate experimental data agreement (if available)

For detailed interpretation, refer to:
- MolProbity: http://molprobity.biochem.duke.edu/
- PDB validation reports: https://www.wwpdb.org/validation/
EOF

echo "✓ Validation summary saved to validation_summary.txt"
cd ..
```

---

## Automated Quality Assessment

### 15. Web-based Validation Services

```python
# Automated submission to web validation services
import requests
import time
import os

def submit_to_molprobity(pdb_file):
    """Submit structure to MolProbity web service"""
    url = "http://molprobity.biochem.duke.edu/index.php"

    # This is a simplified example - actual implementation would need
    # to handle the web form submission properly
    print(f"Manual submission required to: {url}")
    print("Upload your PDB file and download the results")

    return {
        'service': 'MolProbity',
        'url': url,
        'status': 'manual_submission_required'
    }

def submit_to_pdb_redo(pdb_id):
    """Check PDB-REDO for re-refined structure"""
    url = f"https://pdb-redo.eu/db/{pdb_id}"

    try:
        response = requests.get(url)
        if response.status_code == 200:
            return {
                'service': 'PDB-REDO',
                'url': url,
                'status': 'available',
                'message': 'Re-refined structure available'
            }
        else:
            return {
                'service': 'PDB-REDO',
                'status': 'not_available'
            }
    except:
        return {
            'service': 'PDB-REDO',
            'status': 'error'
        }

def check_validation_services(structure_file, pdb_id=None):
    """Check multiple validation services"""
    results = []

    # MolProbity
    results.append(submit_to_molprobity(structure_file))

    # PDB-REDO (if PDB ID available)
    if pdb_id:
        results.append(submit_to_pdb_redo(pdb_id))

    return results
```

### 16. Integration with PyMOL for Visualization

```python
# PyMOL script for structure quality visualization
def generate_pymol_quality_script(structure_file, output_script="quality_check.pml"):
    """Generate PyMOL script for quality visualization"""

    script_content = f"""
# PyMOL Structure Quality Visualization Script
# Generated for: {structure_file}

# Load structure
load {structure_file}, protein

# Basic representation
hide everything
show cartoon, protein
color gray80, protein

# Highlight quality issues
# 1. B-factor coloring
spectrum b, blue_white_red, protein
show sticks, b > 50  # High B-factor regions

# 2. Show missing density (low occupancy)
select low_occ, q < 0.5
show spheres, low_occ
color red, low_occ

# 3. Ramachandran outliers (simplified - would need actual calculation)
# This is a placeholder - real implementation would calculate phi/psi angles
select ramachandran_outliers, none  # Would be populated with actual outliers
color magenta, ramachandran_outliers

# 4. Potential clashes (simplified)
select clashes, none  # Would be populated with actual clashes
show spheres, clashes
color yellow, clashes

# 5. Missing atoms/residues
# Show as gaps in cartoon representation (automatic in PyMOL)

# Create views
set_view (\
    1.000000000,    0.000000000,    0.000000000,\
    0.000000000,    1.000000000,    0.000000000,\
    0.000000000,    0.000000000,    1.000000000,\
    0.000000000,    0.000000000, -100.000000000,\
    0.000000000,    0.000000000,    0.000000000,\
  -20.000000000,  220.000000000,  -20.000000000 )

# Save session
save quality_analysis.pse

# Generate images
png b_factor_view.png, dpi=300
ray 1200, 900
png high_quality_view.png, dpi=300

# Print summary
print "Quality visualization complete!"
print "Red spheres: Low occupancy atoms"
print "Yellow spheres: Potential clashes"
print "Magenta: Ramachandran outliers"
print "Color spectrum: B-factor (blue=low, red=high)"

# Legend
pseudoatom legend1, pos=[10,10,10], color=blue, label="Low B-factor"
pseudoatom legend2, pos=[10,12,10], color=red, label="High B-factor"
pseudoatom legend3, pos=[10,14,10], color=yellow, label="Clashes"
pseudoatom legend4, pos=[10,16,10], color=magenta, label="Rama outliers"

show spheres, legend*
set sphere_scale, 0.5, legend*
"""

    with open(output_script, 'w') as f:
        f.write(script_content)

    print(f"PyMOL script saved to: {output_script}")
    print("Run with: pymol quality_check.pml")

    return output_script
```

---

## Best Practices and Interpretation

### 17. Quality Thresholds and Interpretation Guidelines

```python
class QualityThresholds:
    """Standard quality thresholds for structure assessment"""

    RESOLUTION = {
        'excellent': {'xray': 2.0, 'cryo_em': 3.0},
        'good': {'xray': 2.5, 'cryo_em': 4.0},
        'acceptable': {'xray': 3.5, 'cryo_em': 6.0}
    }

    R_FACTORS = {
        'excellent': {'r_free': 0.20, 'gap': 0.03},
        'good': {'r_free': 0.25, 'gap': 0.05},
        'questionable': {'r_free': 0.30, 'gap': 0.07}
    }

    RAMACHANDRAN = {
        'excellent': {'favored': 98, 'outliers': 0.5},
        'good': {'favored': 95, 'outliers': 1.0},
        'acceptable': {'favored': 90, 'outliers': 2.0}
    }

    GEOMETRY = {
        'bond_rmsd': {'excellent': 0.02, 'good': 0.03, 'poor': 0.05},
        'angle_rmsd': {'excellent': 2.0, 'good': 3.0, 'poor': 5.0}
    }

    B_FACTORS = {
        'wilson_b': {'good': 30, 'acceptable': 50, 'poor': 80},
        'mean_b': {'good': 25, 'acceptable': 40, 'poor': 60}
    }

def interpret_quality_metrics(results):
    """Interpret quality assessment results"""

    interpretation = {
        'overall_grade': 'A',  # Will be determined
        'issues': [],
        'strengths': [],
        'recommendations': []
    }

    # Resolution assessment
    resolution = results.get('metadata', {}).get('resolution')
    method = results.get('metadata', {}).get('method', '').upper()

    if resolution:
        if 'X-RAY' in method:
            if resolution <= QualityThresholds.RESOLUTION['excellent']['xray']:
                interpretation['strengths'].append("Excellent resolution for X-ray structure")
            elif resolution <= QualityThresholds.RESOLUTION['good']['xray']:
                interpretation['strengths'].append("Good resolution")
            else:
                interpretation['issues'].append("Low resolution may limit detailed analysis")

        elif 'ELECTRON MICROSCOPY' in method:
            if resolution <= QualityThresholds.RESOLUTION['excellent']['cryo_em']:
                interpretation['strengths'].append("Excellent resolution for cryo-EM structure")
            elif resolution <= QualityThresholds.RESOLUTION['good']['cryo_em']:
                interpretation['strengths'].append("Good resolution for cryo-EM")
            else:
                interpretation['issues'].append("Resolution limits atomic detail")

    # R-factor assessment
    r_free = results.get('metadata', {}).get('r_free')
    r_work = results.get('metadata', {}).get('r_work')

    if r_free and r_work:
        gap = r_free - r_work

        if r_free <= QualityThresholds.R_FACTORS['excellent']['r_free']:
            interpretation['strengths'].append("Excellent R-factors")
        elif r_free <= QualityThresholds.R_FACTORS['good']['r_free']:
            interpretation['strengths'].append("Good R-factors")
        else:
            interpretation['issues'].append("High R-factors suggest poor model-data agreement")

        if gap > QualityThresholds.R_FACTORS['questionable']['gap']:
            interpretation['issues'].append("Large R-free gap suggests overfitting")

    # Ramachandran assessment
    rama = results.get('ramachandran', {})
    if rama:
        outlier_percent = rama.get('outlier_percent', 0)
        favored_percent = rama.get('favored_percent', 0)

        if (outlier_percent <= QualityThresholds.RAMACHANDRAN['excellent']['outliers'] and
            favored_percent >= QualityThresholds.RAMACHANDRAN['excellent']['favored']):
            interpretation['strengths'].append("Excellent backbone geometry")
        elif outlier_percent <= QualityThresholds.RAMACHANDRAN['acceptable']['outliers']:
            interpretation['strengths'].append("Good backbone geometry")
        else:
            interpretation['issues'].append("Backbone geometry issues detected")
            interpretation['recommendations'].append("Review Ramachandran outliers manually")

    # Completeness assessment
    completeness = results.get('completeness', {})
    if completeness:
        avg_completeness = sum(chain['completeness'] for chain in completeness.values()) / len(completeness)

        if avg_completeness >= 95:
            interpretation['strengths'].append("Model is highly complete")
        elif avg_completeness >= 90:
            interpretation['strengths'].append("Model is reasonably complete")
        else:
            interpretation['issues'].append("Significant missing atoms/residues")
            interpretation['recommendations'].append("Check for missing density in problematic regions")

    # Overall grade determination
    num_issues = len(interpretation['issues'])
    num_strengths = len(interpretation['strengths'])

    if num_issues == 0 and num_strengths >= 3:
        interpretation['overall_grade'] = 'A'
    elif num_issues <= 1 and num_strengths >= 2:
        interpretation['overall_grade'] = 'B'
    elif num_issues <= 2:
        interpretation['overall_grade'] = 'C'
    else:
        interpretation['overall_grade'] = 'D'

    # General recommendations
    if interpretation['overall_grade'] in ['A', 'B']:
        interpretation['recommendations'].append("Structure suitable for most applications")
    elif interpretation['overall_grade'] == 'C':
        interpretation['recommendations'].append("Use with caution for detailed analysis")
        interpretation['recommendations'].append("Consider alternative structures if available")
    else:
        interpretation['recommendations'].append("Structure has significant quality issues")
        interpretation['recommendations'].append("Not recommended for detailed structural analysis")

    return interpretation
```

### 18. Application-Specific Quality Requirements

```python
def assess_fitness_for_purpose(results, application):
    """Assess structure quality for specific applications"""

    applications = {
        'drug_design': {
            'min_resolution': 2.5,
            'max_r_free': 0.25,
            'min_completeness': 95,
            'critical_regions': ['active_site', 'binding_pocket']
        },
        'homology_modeling': {
            'min_resolution': 3.0,
            'max_r_free': 0.30,
            'min_completeness': 90,
            'critical_regions': ['conserved_regions']
        },
        'molecular_dynamics': {
            'min_resolution': 2.8,
            'max_r_free': 0.28,
            'min_completeness': 98,
            'critical_regions': ['all_atoms']
        },
        'functional_studies': {
            'min_resolution': 3.5,
            'max_r_free': 0.35,
            'min_completeness': 85,
            'critical_regions': ['functional_sites']
        }
    }

    if application not in applications:
        return {"error": f"Unknown application: {application}"}

    criteria = applications[application]
    assessment = {"suitable": True, "issues": [], "recommendations": []}

    # Check resolution
    resolution = results.get('metadata', {}).get('resolution')
    if resolution and resolution > criteria['min_resolution']:
        assessment['suitable'] = False
        assessment['issues'].append(f"Resolution ({resolution:.2f} Å) exceeds requirement ({criteria['min_resolution']} Å)")

    # Check R-factors
    r_free = results.get('metadata', {}).get('r_free')
    if r_free and r_free > criteria['max_r_free']:
        assessment['suitable'] = False
        assessment['issues'].append(f"R-free ({r_free:.3f}) exceeds requirement ({criteria['max_r_free']})")

    # Check completeness
    completeness = results.get('completeness', {})
    if completeness:
        avg_completeness = sum(chain['completeness'] for chain in completeness.values()) / len(completeness)
        if avg_completeness < criteria['min_completeness']:
            assessment['suitable'] = False
            assessment['issues'].append(f"Completeness ({avg_completeness:.1f}%) below requirement ({criteria['min_completeness']}%)")

    # Application-specific recommendations
    if application == 'drug_design':
        assessment['recommendations'].extend([
            "Verify binding site geometry and electron density",
            "Check for alternate conformations in binding region",
            "Validate ligand binding poses if co-crystal available"
        ])
    elif application == 'molecular_dynamics':
        assessment['recommendations'].extend([
            "Ensure all atoms are present (especially hydrogens for protonation)",
            "Check for unrealistic conformations that may cause instability",
            "Validate disulfide bonds and metal coordination"
        ])

    return assessment
```

### 19. Comprehensive Validation Report Template

```python
def generate_comprehensive_report(structure_file, results, interpretation, fitness_assessments=None):
    """Generate a comprehensive validation report"""

    report = f"""
COMPREHENSIVE PROTEIN STRUCTURE QUALITY ASSESSMENT
================================================

Structure Information
-------------------
File: {structure_file}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Overall Grade: {interpretation['overall_grade']}

Experimental Data Summary
-----------------------
Method: {results.get('metadata', {}).get('method', 'Unknown')}
Resolution: {results.get('metadata', {}).get('resolution', 'N/A')} Å
R-work: {results.get('metadata', {}).get('r_work', 'N/A')}
R-free: {results.get('metadata', {}).get('r_free', 'N/A')}
Deposition Date: {results.get('metadata', {}).get('deposition_date', 'N/A')}

Model Quality Metrics
-------------------
Ramachandran Analysis:
  • Favored regions: {results.get('ramachandran', {}).get('favored_percent', 'N/A'):.1f}%
  • Allowed regions: {results.get('ramachandran', {}).get('allowed_percent', 'N/A'):.1f}%
  • Outliers: {results.get('ramachandran', {}).get('outlier_percent', 'N/A'):.1f}%

Geometric Validation:
  • Clashes detected: {len(results.get('clashes', []))}
  • Bond length RMSD: {results.get('geometry', {}).get('bond_rmsd', 'N/A')} Å
  • Bond angle RMSD: {results.get('geometry', {}).get('angle_rmsd', 'N/A')}°

Model Completeness:
"""

    # Add completeness details
    for chain_id, comp in results.get('completeness', {}).items():
        report += f"  • Chain {chain_id}: {comp['completeness']:.1f}% complete\n"
        if comp['missing_residues']:
            report += f"    Missing residues: {', '.join(map(str, comp['missing_residues'][:10]))}\n"

    report += f"""
B-factor Analysis:
  • Mean B-factor: {results.get('b_factors', {}).get('mean_b', 'N/A'):.1f} Ų
  • B-factor range: {results.get('b_factors', {}).get('min_b', 'N/A'):.1f} - {results.get('b_factors', {}).get('max_b', 'N/A'):.1f} Ų
  • Standard deviation: {results.get('b_factors', {}).get('std_b', 'N/A'):.1f} Ų

Quality Assessment Summary
------------------------
Strengths:
"""

    for strength in interpretation['strengths']:
        report += f"  ✓ {strength}\n"

    if interpretation['issues']:
        report += "\nIssues Identified:\n"
        for issue in interpretation['issues']:
            report += f"  ⚠ {issue}\n"

    report += "\nRecommendations:\n"
    for rec in interpretation['recommendations']:
        report += f"  • {rec}\n"

    # Add fitness-for-purpose assessments if provided
    if fitness_assessments:
        report += "\nFitness for Purpose Assessment\n"
        report += "-----------------------------\n"
        for application, assessment in fitness_assessments.items():
            suitability = "✓ SUITABLE" if assessment['suitable'] else "✗ NOT SUITABLE"
            report += f"{application.replace('_', ' ').title()}: {suitability}\n"
            if assessment['issues']:
                for issue in assessment['issues']:
                    report += f"  - {issue}\n"

    report += f"""
Quality Score: {results.get('quality_score', 'N/A')}/100

Validation Resources
------------------
For additional validation, consider:
• MolProbity: http://molprobity.biochem.duke.edu/
• PDB Validation Reports: https://www.wwpdb.org/validation/
• SAVES: https://saves.mbi.ucla.edu/
• What Check: https://swift.cmbi.umcn.nl/whatif/

Report Generated by Protein Structure Quality Assessment Pipeline
"""

    return report
```

### 20. Quality Control Checklist

```markdown
## Structure Quality Control Checklist

### Pre-Analysis Setup

- [ ] Verify file format (PDB/mmCIF) and integrity
- [ ] Check for presence of experimental data files (MTZ, MAP)
- [ ] Confirm structure method (X-ray, NMR, Cryo-EM)
- [ ] Identify any special features (ligands, metals, modifications)

### Experimental Data Quality

- [ ] Resolution appropriate for intended use
- [ ] R-factors within acceptable ranges
- [ ] R-free gap < 7%
- [ ] Wilson B-factor reasonable
- [ ] Data completeness > 95% (for X-ray)

### Geometric Validation

- [ ] Ramachandran plot: >95% favored, <2% outliers
- [ ] Bond lengths within 3σ of ideal values
- [ ] Bond angles within 3σ of ideal values
- [ ] No serious atomic clashes (>0.4 Å overlap)
- [ ] Proper stereochemistry for chiral centers

### Model Completeness

- [ ] All expected residues present or justified absences
- [ ] Side chains built where electron density supports
- [ ] Waters placed in reasonable density
- [ ] Metal coordination geometries correct
- [ ] Ligands properly modeled and validated

### B-factor Analysis

- [ ] Reasonable B-factor distribution
- [ ] No atoms with zero occupancy (unless justified)
- [ ] Alternate conformations properly modeled
- [ ] B-factors correlate with structural features

### Final Validation

- [ ] Cross-validation with independent methods
- [ ] Comparison with homologous structures
- [ ] Biological relevance assessment
- [ ] Documentation of limitations and uncertainties

### Application-Specific Checks

#### For Drug Design:

- [ ] Active site well-defined and complete
- [ ] Binding pocket accessible and properly shaped
- [ ] Cofactor/metal sites correctly coordinated

#### For Molecular Dynamics:

- [ ] All atoms present (including hydrogens)
- [ ] No steric clashes that would cause instability
- [ ] Proper protonation states assigned

#### For Functional Studies:

- [ ] Functionally important regions well-resolved
- [ ] Conformational states biologically relevant
- [ ] Quaternary structure validated if applicable
```

This comprehensive guide provides you with all the essential tools and knowledge needed to
thoroughly assess protein structure quality. The combination of automated scripts, manual validation
techniques, and interpretation guidelines will help ensure you're working with reliable structural
data for your research applications.

Key takeaways:

1. **Always consider the intended application** - different uses have different quality requirements
2. **Use multiple validation methods** - no single metric tells the complete story
3. **Pay attention to experimental details** - resolution, R-factors, and data collection methods
   matter
4. **Validate geometry rigorously** - poor geometry indicates modeling problems
5. **Document limitations clearly** - be transparent about quality issues and their implications

The tools and scripts provided can be adapted for your specific needs and integrated into automated
pipelines for high-throughput analysis.
