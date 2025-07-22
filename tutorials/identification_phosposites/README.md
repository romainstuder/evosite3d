# OpenMS Tutorial: Phosphosite Identification in Human Samples

## Theoretical concept

Phosphorylation is a fundamental post-translational modification where phosphate groups are
covalently attached to proteins, typically on serine (S), threonine (T), and tyrosine (Y) residues.
This reversible modification serves as a critical regulatory mechanism controlling protein function,
localisation, stability, and interactions across virtually all cellular processes including signal
transduction, metabolism, and cell cycle progression. The detection and quantification of
phosphorylation events has been revolutionised by proteomics approaches utilising mass spectrometry
(MS), which enable genome-wide identification and characterisation of phosphorylation sites. Modern
phosphoproteomic workflows combine enrichment strategies such as immobilised metal affinity
chromatography (IMAC) or titanium dioxide (TiO2) to selectively capture phosphopeptides, followed by
high-resolution mass spectrometry analysis that can precisely localise phosphorylation sites and
quantify their stoichiometry across different biological conditions. These technological advances
have transformed our understanding of phosphorylation networks and their dysregulation in diseases,
making phosphoproteomics an indispensable tool in systems biology and biomedical research.

## Overview

This tutorial demonstrates how to identify phosphorylation sites (phosphosites) in human proteins
using OpenMS with publicly available mass spectrometry data. We will use a complete workflow with
OpenMS from raw data processing to phosphosite localisation.

https://openms.de/

```
OpenMS offers an open-source C++ library (+ Python bindings) for LC/MS data management, analysis and
visualization. It empowers rapid development of mass spectrometry related software. OpenMS is freely
available under the three clause BSD license and runs under Windows, macOS and Linux. The OpenMS
members have a strong commitment to creating an open, inclusive, and positive community. Please read
the OpenMS Code of Conduct for guidance on how to interact with others in a way that makes the
community thrive.
```

## Prerequisites

### Software Requirements

- OpenMS 3.0+ (with TOPP tools)

Optional:

- TOPPView (for visualization)
- SearchGUI and PeptideShaker (optional, for alternative search)

### Installation

#### macOS

```shell

# Method 1: Using Homebrew
brew install openms
# Important, you need to add the $DYLD_LIBRARY_PATH path:
export DYLD_LIBRARY_PATH="/opt/homebrew/opt/openms/lib:$DYLD_LIBRARY_PATH"


# Method 2: Download pre-built binaries
# Visit: https://openms.de/downloads/
# Download the macOS installer (.dmg file)

# Method 3: Build from source using conda
conda create -n openms-env python=3.9
conda activate openms-env
conda install -c bioconda openms
```

#### Ubuntu/Debian

```shell
sudo apt-get update
sudo apt-get install openms openms-doc
```

#### Windows

```shell
# Download installer from: https://openms.de/downloads/
# Or use conda:
conda create -n openms-env python=3.9
conda activate openms-env
conda install -c bioconda openms
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

## Step 2: File Conversion and Preprocessing

### Step 2.1: Define files used downstream

```shell
SPECIES="yeast"

case "$SPECIES" in
  "human")
    FILES=(
        "20120302_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_EGF5_rep1_pH6"
        "20111219_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_Control_rep1_pH6"
    )
    ;;
  "yeast")
    FILES=(
        "34160_e02555_ML040_phos_StressPool_PFP3"
    )
    ;;
esac
```

### Step 2.2: Convert Raw Files to mzML

First, we need to convert the raw dataset from MS to another format, mzML.

Raw files (.raw) are proprietary binary formats produced by mass spectrometers that contain the
original instrument data including scan metadata, ion intensities, and timing information. The mzML
format is an open, standardized XML-based format that stores the same mass spectrometry data in a
vendor-neutral way, making it compatible with various analysis software tools. Converting from raw
to mzML ensures data accessibility and reproducibility across different platforms while preserving
all essential spectral information.

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

### Finally, we can inspect the files:

```shell
for base_name in "${FILES[@]}"; do
    mzml_file="${base_name}.mzML"
    FileConverter -in "$mzml_file"
done
```

## Step 3: Perform peak picking (batch processing)

Peak picking is the process of identifying and extracting meaningful signals (peaks) from continuous
mass spectrometry data by distinguishing real analyte signals from background noise and baseline
fluctuations. This step converts profile-mode data into centroided data by determining peak centers,
intensities, and boundaries using signal-to-noise ratio calculations. Proper peak picking is crucial
for accurate mass measurements and downstream identification algorithms, as it directly affects the
quality of peptide spectrum matches.

https://openms.de/current_doxygen/html/TOPP_PeakPickerHiRes.html
https://openms.readthedocs.io/en/latest/getting-started/types-of-topp-tools/picking-peaks.html

```shell
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

## Step 4: Database Search

### Step 4.1: Prepare Protein Database

```shell
SPECIES="yeast"

case "$SPECIES" in
  "human")
    TAXID="9606"
    PROTEOME="UP000005640"
    ;;
  "yeast")
    TAXID="559292"
    PROTEOME="UP000002311"
    ;;
esac

# Download human proteome from UniProt
wget "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/${PROTEOME}/${PROTEOME}_${TAXID}.fasta.gz"

# Decompress
gunzip "${PROTEOME}_${TAXID}.fasta.gz"

# Rename for clarity
mv "${PROTEOME}_${TAXID}.fasta" "${SPECIES}_proteome.fasta"
```

### Step 4.3: Create Decoy Database

Decoy databases contain reversed or shuffled protein sequences that serve as negative controls for
statistical validation of peptide identifications in database searches. When searching against a
combined target-decoy database, false matches to decoy sequences estimate the false discovery rate
(FDR) of identifications. This approach allows researchers to control identification confidence by
setting FDR thresholds, typically 1-5%, ensuring that the majority of reported identifications are
genuine matches rather than random coincidences.

https://www.openms.org/documentation/html/TOPP_DecoyDatabase.html

```shell
DecoyDatabase -in "${SPECIES}_proteome.fasta" \
    -out "${SPECIES}_proteome_decoy.fasta"" \
    -decoy_string "DECOY_" \
    -method reverse
```

#TODO: Warning: Only one FASTA input file was provided, which might not contain contaminants. You
probably want to have them! Just add the contaminant file to the input file list 'in'.

### Step 4.2: Configure Search Parameters

Database searching involves comparing experimental peptide mass spectra against theoretical spectra
generated from a protein sequence database to identify which peptides produced the observed
fragmentation patterns. Search engines like MS-GF+ and COMET score the similarity between
experimental and theoretical spectra, considering factors like precursor mass accuracy, fragment ion
matches, and modification states. This computational approach enables high-throughput identification
of thousands of peptides from complex proteome samples.

We will use the following TOPP INI file (in XML format):
[search_params.xml](./data/search_params.xml)

Ref: https://openms.readthedocs.io/en/latest/getting-started/topp-tools.html

### Step 4.3: Search for phosphopeptides, using batch process with MS-GF+

https://www.openms.org/documentation/html/TOPP_MSGFPlusAdapter.html

```shell
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

```shell
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

False Discovery Rate (FDR) control in phosphoproteomics ensures that a specified percentage
(typically 1%) of reported phosphopeptide identifications are expected to be incorrect matches. FDR
calculation uses the ratio of decoy hits to target hits, multiplied by 100, to estimate the
proportion of false positives in the identification list. Controlling FDR at the peptide-spectrum
match (PSM) level is particularly important for phosphoproteomics due to the added complexity of
modification site assignment and the lower abundance of phosphorylated species.

https://openms.de/documentation/html/TOPP_FalseDiscoveryRate.html

```shell
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

Identification filtering applies quality criteria to remove low-confidence peptide matches based on
multiple parameters including identification scores, peptide length, precursor charge state, and
mass accuracy. These filters eliminate identifications that are likely to be incorrect due to poor
spectral quality, unusual peptide properties, or low search engine confidence scores. Proper
filtering significantly improves the reliability of downstream analyses by retaining only
high-quality identifications that meet stringent analytical criteria.

```shell
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

Peptide indexing maps identified peptides back to their parent proteins in the database,
establishing the relationship between peptide sequences and protein accession numbers. This process
handles cases where peptides match multiple proteins (shared peptides) and ensures proper protein
grouping for quantification and functional analysis. Indexing also validates that identified
peptides conform to expected enzyme cleavage patterns and protein sequence boundaries, providing
additional confidence in the identification results.

```shell
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

Phosphosite localisation determines the exact amino acid residue(s) within a peptide that carry the
phosphorylation modification, which is critical since multiple serine, threonine, or tyrosine
residues may be present in the same peptide. Algorithms like PhosphoRS or Ascore analyse
fragmentation patterns to calculate localisation probabilities for each potential modification site.
High-confidence site localisation (typically >75% probability) is essential for understanding the
biological significance of specific phosphorylation events and their regulatory roles.

https://openms.de/documentation/TOPP_PhosphoScoring.html

```shell
for base_name in "${FILES[@]}"; do
    input_file="${base_name}_picked.mzML"
    index_file="${base_name}_indexed.idXML"
    output_file="${base_name}_phosphoscoring_results.idXML"
    echo "Processing $input_file..."

    # Use PhosphoScoring for site localisation
    PhosphoScoring -in "$input_file" \
        -id "$index_file" \
        -out "$output_file" \
        -threads 4

    echo "Completed: $output_file"
done
```

## Step 7: Results Export and Analysis

### Step 7.1: Export to Text Format

```shell
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

```shell
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

### Visualisation with TOPPView

```shell
# Launch TOPPView for visual inspection
TOPPView picked.mzML localized.idXML
```

### Quality Metrics

```shell
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

```shell
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

```shell
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
2. **localisation confidence**: 60-80% with score >0.75
3. **Common modifications**: pSer (85%), pThr (12%), pTyr (3%)

### Result Interpretation

- **localisation score >0.75**: High confidence
- **localisation score 0.5-0.75**: Medium confidence
- **localisation score <0.5**: Low confidence (ambiguous)

## Troubleshooting

### Common Issues

1. **Low identification rate**: Check search parameters, enzyme specificity
2. **Poor localisation**: Verify fragment mass tolerance, neutral losses
3. **High FDR**: Increase database size, check decoy generation

### Performance Tips

- Use SSD storage for faster I/O
- Allocate sufficient RAM (8GB+ recommended)
- Consider parallel processing for multiple files

## Advanced Applications

### Quantitative Phosphoproteomics

```shell
# For TMT/iTRAQ quantification
IsobaricAnalyzer -in picked.mzML \
    -out quantified.idXML \
    -extraction:select_activation HCD \
    -extraction:reporter_mass_shift 0.002
```

### Time-Course Analysis

```shell
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
