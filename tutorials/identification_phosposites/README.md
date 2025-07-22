# OpenMS Tutorial: Phosphosite Identification in Human Samples

## Overview

This tutorial demonstrates how to identify phosphorylation sites (phosphosites) in human proteins
using OpenMS with publicly available mass spectrometry data. We'll use a complete workflow from raw
data processing to phosphosite localization.

## Prerequisites

### Software Requirements

- OpenMS 3.0+ (with TOPP tools)
- TOPPView (for visualization)
- SearchGUI and PeptideShaker (optional, for alternative search)

### Installation

#### macOS

```bash
# Method 1: Download pre-built binaries
# Visit: https://openms.de/downloads/
# Download the macOS installer (.dmg file)

# Method 3: Build from source using conda
conda create -n openms-env python=3.9
conda activate openms-env
conda install -c bioconda openms

# Method 3: Using Homebrew (if available)
brew install openms
```

#### Ubuntu/Debian

```bash
sudo apt-get update
sudo apt-get install openms openms-doc
```

#### Windows

```bash
# Download installer from: https://openms.de/downloads/
# Or use conda:
conda create -n openms-env python=3.9
conda activate openms-env
conda install -c bioconda openms
```

```shell
export DYLD_LIBRARY_PATH="/opt/homebrew/opt/openms/lib:$DYLD_LIBRARY_PATH"
```

### Required Files

- Human protein database (UniProt)
- Phosphorylation modification database
- Public phosphoproteomics dataset

## Step 1: Data Acquisition

### Download Public Dataset - PXD000612

We'll use the "Ultra-deep human phosphoproteome" dataset from HeLa cells:

**Dataset Information:**

- **PXD ID**: PXD000612
- **Title**: Ultra-deep human phosphoproteome reveals different regulatory nature of Tyr and
  Ser/Thr-based signaling
- **Organism**: Homo sapiens (HeLa cells)
- **Instrument**: Q Exactive
- **Contains**: >50,000 phosphorylated peptides

**Recommended Files to Download:**

For a tutorial, start with these representative files (smaller subset):

```
# Access the dataset at: https://www.ebi.ac.uk/pride/archive/projects/PXD000612

# Recommended starting files (choose 2-3 for tutorial):
# 1. Control/baseline samples:
#    - HeLa_Control_1.raw
#    - HeLa_Control_2.raw

# 2. Stimulated samples (if doing comparative analysis):
#    - HeLa_EGF_1.raw
#    - HeLa_EGF_2.raw

# 3. Alternative: Pick any 2-3 .raw files from the file list
#    (The exact filenames may vary - check the FTP download section)

# Additional files that might be helpful:
# - mqpar.xml (MaxQuant parameters used in original study)
# - Any .fasta database files provided
# - Peak lists (.mgf files) if available
```

**File Selection Tips:**

- Start with 2-3 files to test the workflow (each file is several GB)
- Look for files with clear naming (Control, Treatment, replicate numbers)
- Avoid downloading the entire dataset initially (>100 files, several TB)
- Focus on technical replicates from the same condition first

# https://www.ebi.ac.uk/pride/archive/projects/PXD035029

**Download Methods:**

```shell
# Method 1: Direct download from PRIDE web interface
# 1. Go to https://www.ebi.ac.uk/pride/archive/projects/PXD000612
# 2. Click "FTP Download"
# 3. Select representative .raw files

# Method 2: Automated download using wget (recommended)
BASE_URL="ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2014/08/PXD000612"
FILES=(
    "20120302_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_EGF5_rep1_pH6.raw"
    "20120302_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_EGF5_rep2_pH6.raw"
    "20120302_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_EGF5_rep3_pH6.raw"
    "20111219_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_Control_rep1_pH6.raw"
    "20111219_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_Control_rep2_pH6.raw"
    "20111219_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_Control_rep3_pH6.raw"
)
wget -c -t 3 "${FILES[@]/#/$BASE_URL/}"

# Alternative: Download subset for faster tutorial
BASE_URL="ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2014/08/PXD000612"
TUTORIAL_FILES=(
    "20120302_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_EGF5_rep1_pH6.raw"
    "20111219_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_Control_rep1_pH6.raw"
)
wget -c -t 3 "${TUTORIAL_FILES[@]/#/$BASE_URL/}"
```

### Prepare Protein Database

```bash
# Download human proteome from UniProt
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002311/UP000002311_559292.fasta.gz

# Decompress
gunzip UP000005640_9606.fasta.gz
gunzip UP000002311_559292.fasta.gz

# Rename for clarity
mv UP000005640_9606.fasta human_proteome.fasta
mv UP000002311_559292.fasta yeast_proteome.fasta
```

## Step 2: File Conversion and Preprocessing

### Define files used downstream

```shell
FILES=(
    "20120302_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_EGF5_rep1_pH6"
    "20111219_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_Control_rep1_pH6"
)

FILES=(
    "34160_e02555_ML040_phos_StressPool_PFP3"
)

```

### Step 2.1: Convert Raw Files to mzML

```shell
for base_name in "${FILES[@]}"; do
    raw_file="${base_name}.raw"
    mzml_file="${base_name}.mzML"

    echo "Processing $input_file..."

    FileConverter -in "$raw_file" \
       -out "$mzml_file" \
       -RawToMzML:ThermoRaw_executable "../ThermoRawFileParser1.4.5/ThermoRawFileParser.exe"
done
```

### We can inspect the files:

```shell
FileInfo -in 20111219_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_Control_rep1_pH6.mzML
FileInfo -in 20120302_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_EGF5_rep1_pH6.mzML
FileInfo -in 34160_e02555_ML040_phos_StressPool_PFP3.mzML
```

### Step 2.2: Perform peak picking (batch processing)

https://openms.de/current_doxygen/html/TOPP_PeakPickerHiRes.html

```bash
for base_name in "${FILES[@]}"; do
    mzml_file="${base_name}.mzML"
    picked_file="${base_name}_picked.mzML"

    echo "Processing mzml_file..."

    PeakPickerHiRes -in "$mzml_file" \
       -out "$picked_file" \
       -algorithm:signal_to_noise 1.0 \
       -threads 4
done
```

## Step 3: Database Search Configuration

### Step 3.1: Create Decoy Database

https://www.openms.org/documentation/html/TOPP_DecoyDatabase.html

```bash
DecoyDatabase -in human_proteome.fasta \
    -out human_proteome_decoy.fasta \
    -decoy_string "DECOY_" \
    -method reverse
```

```shell
DecoyDatabase -in yeast_proteome.fasta \
    -out yeast_proteome_decoy.fasta \
    -decoy_string "DECOY_" \
    -method reverse
```

#TODO: Warning: Only one FASTA input file was provided, which might not contain contaminants. You
probably want to have them! Just add the contaminant file to the input file list 'in'.

### Step 3.2: Configure Search Parameters

Create a search parameter file (`search_params.xml`):

```xml

<PARAMETERS>
<NODE name="MSGFPlusAdapter">
  <NODE name="1">
    <!-- Use MS-GF+ for phospho searches -->
    <ITEM name="instrument" value="Q_Exactive" type="string"/>
    <ITEM name="enzyme" value="Trypsin/P" type="string"/>
    <ITEM name="max_missed_cleavages" value="2" type="int"/>
    <ITEM name="precursor_mass_tolerance" value="20.0" type="float"/>
    <ITEM name="precursor_error_units" value="ppm" type="string"/>
    <ITEM name="max_mods" value="2" type="int"/>
    <!-- Fixed modifications -->
    <ITEM name="fixed_modifications" value="Carbamidomethyl (C)" type="mods"/>
    <!-- Variable modifications for phosphorylation -->
    <ITEM name="variable_modifications" value="Oxidation (M), Phospho (STY)" type="mods"/>
  </NODE>
</NODE>
</PARAMETERS>
```

## Step 4: Database Search

### Step 4.1: Search for phosphopeptides, using batch process with MS-GF+

https://www.openms.org/documentation/html/TOPP_MSGFPlusAdapter.html

```bash
SPECIES="yeast"
for base_name in "${FILES[@]}"; do
    input_file="${base_name}_picked.mzML"
    output_file="${base_name}_search_results.idXML"

    echo "Processing $input_file..."

    MSGFPlusAdapter -in "$input_file" \
        -out "$output_file" \
        -database "./${SPECIES}_proteome_decoy.fasta" \
        -executable ../MSGFPlus/MSGFPlus.jar \
        -ini search_params.xml \
        -java_memory 8000

    echo "Completed: $output_file"
done
```

Alternative with COMET:

```bash
SPECIES="yeast"
for base_name in "${FILES[@]}"; do
    input_file="${base_name}_picked.mzML"
    output_file="${base_name}_search_results.idXML"

    echo "Processing $input_file..."

    CometAdapter -in "$input_file" \
        -out "$output_file" \
        -database "./${SPECIES}_proteome_decoy.fasta" \
        -comet_executable comet.exe \
        -ini search_params.xml

    echo "Completed: $output_file"
done
```

## Step 5: Post-Processing and Filtering

### Step 5.1: Apply False Discovery Rate control

```bash
for base_name in "${FILES[@]}"; do
    input_file="${base_name}_search_results.idXML"
    fdr_file="${base_name}_fdr.idXML"

    echo "Processing $input_file..."

    FalseDiscoveryRate -in "$input_file" \
    -out "$fdr_file" \
    -PSM true \
    -protein false \
    -FDR:PSM 0.01

    echo "Completed: $output_file"
done
```

### Step 5.2: Filter identifications by score and quality

```bash
for base_name in "${FILES[@]}"; do
    fdr_file="${base_name}_fdr.idXML"
    filtered_file="${base_name}_fdr_filtered.idXML"

    echo "Filtering identifications..."

    # Filter identifications by score and quality
    IDFilter -in "$fdr_file" \
        -out "$filtered_file" \
        -score:peptide 0.01 \
        -remove_decoys \
        -precursor:length '8:50' \
        -precursor:charge '2:5'

    echo "Completed: $output_file"
done
```

### ### Step 5.3: Peptide Indexing

```bash
for base_name in "${FILES[@]}"; do
    filtered_file="${base_name}_fdr_filtered.idXML"
    indexed_file="${base_name}_indexed.idXML"

    echo "Processing $input_file..."

    # Index peptides to proteins
    PeptideIndexer -in "$filtered_file" \
        -out "$indexed_file" \
        -fasta "./${SPECIES}_proteome_decoy.fasta" \
        -enzyme:specificity full

    echo "Completed: $output_file"
done
```

## Step 6: Phosphosite Localisation

https://openms.de/documentation/TOPP_PhosphoScoring.html

```bash
for base_name in "${FILES[@]}"; do
    input_file="${base_name}_picked.mzML"
    index_file="${base_name}_indexed.idXML"
    output_file="${base_name}_phosphoscoring_results.idXML"
    echo "Processing $input_file..."

    # Use PhosphoScoring for site localization
    PhosphoScoring -in "$input_file" \
        -id "$index_file" \
        -out "$output_file" \
        -threads 4

    echo "Completed: $output_file"
done
```

## Step 7: Results Export and Analysis

### Step 7.1: Export to Text Format

```bash
# Search for phosphopeptides, using batch process with MS-GF+
for base_name in "${FILES[@]}"; do
    input_file="${base_name}_phosphoscoring_results.idXML"
    output_file="${base_name}_phosphosite_results.tsv"
    echo "Processing $input_file..."

    # Export identification results
    TextExporter -in "$input_file" \
        -out "$output_file" \
        -out_type tsv \
        -id:peptides_only \
        -id:first_dim_rt
    echo "Completed: $output_file"
done
```

# More stringent filtering for phospho analysis

```
IDFilter -in fdr_filtered.idXML \
    -out filtered.idXML \
    -score:pep 0.01 \
    -remove_decoys \
    -precursor:length '8:50' \
    -precursor:charge '2:5' \
    -whitelist:modifications 'Phospho (STY)' \  # Only phosphorylated peptides
    -score:prot 0.05                            # Protein-level score threshold
```

### Step 7.2: Export Phosphosite Summary

```bash
# Map peptide identifications to features for quantification
for base_name in "${FILES[@]}"; do
    id_file="${base_name}_indexed.idXML"  # Use indexed results
    feature_file="${base_name}_features.featureXML"
    output_file="${base_name}_quantified.featureXML"

    IDMapper -id "$id_file" \
        -in "$feature_file" \
        -out "$output_file" \
        -rt_tolerance 30.0 \
        -mz_tolerance 0.01
done
```

## Step 8: Quality Control and Validation

### Visualization with TOPPView

```bash
# Launch TOPPView for visual inspection
TOPPView picked.mzML localized.idXML
```

### Quality Metrics

```bash
# Generate QC report
for base_name in "${FILES[@]}"; do
    input_file="${base_name}_picked.mzML"
    qc_report="${base_name}_qc_report.qcML"
    phosphoscoring_results="${base_name}_phosphoscoring_results.idXML"

    QCCalculator -in "$input_file" \
        -out "$qc_report" \
        -id "$phosphoscoring_results"
done
```

## Step 9: Quantitative Comparison (Control vs EGF5)

### Step 9.1: Label-Free Quantification

```bash
# After obtaining all search results, perform feature detection for quantification
for base_name in "${FILES[@]}"; do
    input_file="${base_name}_picked.mzML"
    output_file="${base_name}_features.featureXML"

    echo "Feature detection for $input_file..."
    FeatureFinderCentroided -in "$input_file" \
        -out "$output_file" \
        -algorithm:mass_trace:mz_tolerance 0.005 \
        -algorithm:mass_trace:min_spectra 5 \
        -algorithm:isotopic_pattern:mz_tolerance 0.005
done
```

### Link Features to Identifications

```bash
# Map peptide identifications to features for quantification
for base_name in "${FILES[@]}"; do
    id_file="${base_name}_phosphoscoring_results.idXML"
    feature_file="${base_name}_features.featureXML"
    output_file="${base_name}_quantified.featureXML"

    IDMapper -id "$id_file" \
        -in "$feature_file" \
        -out "$output_file" \
        -rt_tolerance 30.0 \
        -mz_tolerance 0.01
done
```

## Expected Results

### Typical Outputs

1. **Phosphosite identifications**: 1000-5000 sites (dataset dependent)
2. **Localization confidence**: 60-80% with score >0.75
3. **Common modifications**: pSer (85%), pThr (12%), pTyr (3%)

### Result Interpretation

- **Localization score >0.75**: High confidence
- **Localization score 0.5-0.75**: Medium confidence
- **Localization score <0.5**: Low confidence (ambiguous)

## Troubleshooting

### Common Issues

1. **Low identification rate**: Check search parameters, enzyme specificity
2. **Poor localization**: Verify fragment mass tolerance, neutral losses
3. **High FDR**: Increase database size, check decoy generation

### Performance Tips

- Use SSD storage for faster I/O
- Allocate sufficient RAM (8GB+ recommended)
- Consider parallel processing for multiple files

## Advanced Applications

### Quantitative Phosphoproteomics

```bash
# For TMT/iTRAQ quantification
IsobaricAnalyzer -in picked.mzML \
    -out quantified.idXML \
    -extraction:select_activation HCD \
    -extraction:reporter_mass_shift 0.002
```

### Time-Course Analysis

```bash
# Process multiple time points
for timepoint in T0 T15 T30 T60; do
    # Process each timepoint
    echo "Processing $timepoint"
    # Run complete pipeline
done
```

## References and Further Reading

1. OpenMS Documentation: https://openms.readthedocs.io/
2. PRIDE Database: https://www.ebi.ac.uk/pride/
3. Phosphorylation databases: PhosphoSitePlus, PhosphoELM
4. Visualization tools: TOPPView, Skyline

## Conclusion

This tutorial provides a complete workflow for phosphosite identification using OpenMS. The pipeline
can be adapted for different experimental designs and extended with additional analysis steps.
Remember to validate results with orthogonal methods and consider biological context in
interpretation.
