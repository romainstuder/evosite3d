# OpenMS Tutorial: Phosphosite Identification in Yeast Samples

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

This tutorial demonstrates how to identify phosphorylation sites (phosphosites) in yeast proteins
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

The tutorial will follow eight different steps:

1. [Data Acquisition](#dataaquisition)
2. [File Conversion](#fileconversion)
3. [Peak picking](#peakpicking)
4. [Database Search](#databasesearch)
5. [Post-Processing and Filtering](#postprocessing)
6. [Phosphosite Localisation](#phosphositeslocalisation)
7. [Results Export and Analysis](#resultsexport)

## Prerequisites

### Software Requirements

- OpenMS 3.0+ (with TOPP tools)
-

Optional:

- TOPPView (for visualization)
- SearchGUI and PeptideShaker (optional, for alternative search)

### Installation of OpenMS

#### macOS

https://openms.readthedocs.io/en/latest/about/installation/installation-on-macos.html#install-via-macos-installer

Method 1: Download pre-built binaries

```shell
# Visit: https://abibuilder.cs.uni-tuebingen.de/archive/openms/OpenMSInstaller/release/latest/
# Download the macOS installer (.pkg file) for Silicon or x66 and install it.
export PATH="/Applications/OpenMS-3.4.1/bin:$PATH"
```

Method 2: Using Homebrew (still WIP)

```shell
# Method 2: Using Homebrew (WIP)
#brew tap brewsci/bio
#brew install openms
# Important, you need to add the $DYLD_LIBRARY_PATH path:
#export DYLD_LIBRARY_PATH="/opt/homebrew/opt/openms/lib:$DYLD_LIBRARY_PATH"
```

Method 3: Build from source using conda

```shell
conda create -n openms-env python=3.9
conda activate openms-env
conda install -c bioconda openms
```

#### GNU Linux - Ubuntu/Debian

https://openms.readthedocs.io/en/latest/about/installation/installation-on-gnu-linux.html

```shell
# Download installer from: https://abibuilder.cs.uni-tuebingen.de/archive/openms/OpenMSInstaller/release/latest/

# or:
sudo apt-get update
sudo apt-get install openms openms-doc
```

#### Windows

https://openms.readthedocs.io/en/latest/about/installation/installation-on-windows.html

```shell
# Download installer from: https://abibuilder.cs.uni-tuebingen.de/archive/openms/OpenMSInstaller/release/latest/
# Or use conda:
conda create -n openms-env python=3.9
conda activate openms-env
conda install -c bioconda openms
```

### Installation of ThermoRawFileParser (used in FileConverter)

The file `ThermoRawFileParser.exe` can be retrieved from here:
https://github.com/CompOmics/ThermoRawFileParser/releases/ The release used in this tutorial is
`1.4.5` The file needs to be in an available path.

### Installation of MSGFPlusAdapter (used in MSGFPlusAdapter)

The file `MSGFPlus.jar` can be retrieved from here: https://github.com/MSGFPlus/msgfplus/releases
The release used in this tutorial is `2024.03.26` The file needs to be in an available path.

### Installation of Comet (used in CometAdapter)

The file `comet.exe` can be retrieved from here: https://github.com/UWPR/Comet/releases The release
used in this tutorial is `2025.02.0` The file needs to be in an available path.

## Step 1: Data Acquisition <a name="dataaquisition"></a>

### Required Files

- Human protein database (UniProt)
- Phosphorylation modification database
- Public phosphoproteomics dataset
-

### Download Public Dataset - PXD035029

This tutorial demonstrates how to identify phosphorylation sites in Saccharomyces cerevisiae using
mass spectrometry data from the PXD035029 dataset. This dataset represents an ultra-deep yeast
phosphoproteome study that investigated stress-responsive phosphorylation networks across 101
environmental and chemical perturbations. Dataset Information

- Dataset: PXD035029 - "The regulatory landscape of the yeast phosphoproteome"
- Organism: Saccharomyces cerevisiae (Baker's yeast)
- Study: Ultra-deep reference phosphoproteomic DDA data
- Coverage: Phosphosites on 59% of the yeast proteome
- Scope: 18% of the proteome harbors regulated phosphosites within 5 minutes of stress exposure
- Reference Database: SGD (Saccharomyces Genome Database)
- Link: https://www.ebi.ac.uk/pride/archive/projects/PXD035029

**File Selection Tips:**

- Start with 2-3 files to test the workflow (each file is several GB)
- Look for files with clear naming (Control, Treatment, replicate numbers)
- Avoid downloading the entire dataset initially (>100 files, several TB)
- Focus on technical replicates from the same condition first

**Download Methods:**

```shell
# Method 1: Direct download from PRIDE web interface
# 1. Go to https://www.ebi.ac.uk/pride/archive/projects/PXD035029
# 2. Click "FTP Download"
# 3. Select representative .raw files

# Method 2: Automated download using wget (recommended)
BASE_URL="ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2023/08/PXD035029"
TUTORIAL_FILES=(
    "34160_e02555_ML040_phos_StressPool_PFP3.raw"
)
wget -c -t 3 "${TUTORIAL_FILES[@]/#/$BASE_URL/}"
```

## Step 2: File Conversion <a name="fileconversion"></a>

### Step 2.1: Define files used downstream

```shell
SPECIES="yeast"

case "$SPECIES" in
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

Link: https://openms.de/documentation/html/TOPP_FileConverter.html

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

Finally, we can inspect the files:

Link: https://openms.de/documentation/html/TOPP_FileInfo.html

```shell
for base_name in "${FILES[@]}"; do
    mzml_file="${base_name}.mzML"
    FileInfo -in "$mzml_file"
done
```

## Step 3: Peak picking <a name="peakpicking"></a>

Peak picking is the process of identifying and extracting meaningful signals (peaks) from continuous
mass spectrometry data by distinguishing real analyte signals from background noise and baseline
fluctuations. This step converts profile-mode data into centroided data by determining peak centers,
intensities, and boundaries using signal-to-noise ratio calculations. Proper peak picking is crucial
for accurate mass measurements and downstream identification algorithms, as it directly affects the
quality of peptide spectrum matches.

Links: https://openms.de/current_doxygen/html/TOPP_PeakPickerHiRes.html

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

## Step 4: Database Search for Peptide Identification <a name="dataaquisition"></a>

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

# Download proteome from UniProt
wget "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/${PROTEOME}/${PROTEOME}_${TAXID}.fasta.gz"

# Decompress
gunzip "${PROTEOME}_${TAXID}.fasta.gz"

# Rename for clarity
mv "${PROTEOME}_${TAXID}.fasta" "${SPECIES}_proteome.fasta"
```

### Step 4.2: Create Decoy Database

Decoy databases contain reversed or shuffled protein sequences that serve as negative controls for
statistical validation of peptide identifications in database searches. When searching against a
combined target-decoy database, false matches to decoy sequences estimate the false discovery rate
(FDR) of identifications. This approach allows researchers to control identification confidence by
setting FDR thresholds, typically 1-5%, ensuring that the majority of reported identifications are
genuine matches rather than random coincidences.

Link: https://www.openms.org/documentation/html/TOPP_DecoyDatabase.html

```shell
DecoyDatabase -in "${SPECIES}_proteome.fasta" \
    -out "${SPECIES}_proteome_decoy.fasta"" \
    -decoy_string "DECOY_" \
    -method reverse
```

### Step 4.3: Configure Search Search for phosphopeptides, using batch process with MS-GF+

Database searching involves comparing experimental peptide mass spectra against theoretical spectra
generated from a protein sequence database to identify which peptides produced the observed
fragmentation patterns. Search engines like MS-GF+ and COMET score the similarity between
experimental and theoretical spectra, considering factors like precursor mass accuracy, fragment ion
matches, and modification states. This computational approach enables high-throughput identification
of thousands of peptides from complex proteome samples.

We will use the following TOPP INI file (in XML format) to define the parameters for the search:
[search_params.xml](./data/search_params.xml)

Link: https://openms.readthedocs.io/en/latest/getting-started/topp-tools.html
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

Link: https://www.openms.org/documentation/html/TOPP_CometAdapter.html

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

### Step 4.5: Peptide Indexing

Peptide indexing maps identified peptides back to their parent proteins in the database,
establishing the relationship between peptide sequences and protein accession numbers. This process
handles cases where peptides match multiple proteins (shared peptides) and ensures proper protein
grouping for quantification and functional analysis. Indexing also validates that identified
peptides conform to expected enzyme cleavage patterns and protein sequence boundaries, providing
additional confidence in the identification results.

Link: https://openms.de/documentation/html/TOPP_PeptideIndexer.html

```shell
for base_name in "${FILES[@]}"; do
    input_file="${base_name}_search_results.idXML"
    indexed_file="${base_name}_indexed.idXML"

    echo "Processing $input_file..."

    # Index peptides to proteins
    PeptideIndexer -in "${input_file}" \
        -out "${indexed_file}" \
        -fasta "./${SPECIES}_proteome_decoy.fasta" \
        -enzyme:specificity full

    echo "PeptideIndexer Completed for: $output_file"
done
```

## Step 5: Post-processing <a name="postprocessing"></a>

### Step 5.1: Apply False Discovery Rate (FDR) Control

False Discovery Rate (FDR) control in phosphoproteomics ensures that a specified percentage
(typically 1%) of reported phosphopeptide identifications are expected to be incorrect matches. FDR
calculation uses the ratio of decoy hits to target hits, multiplied by 100, to estimate the
proportion of false positives in the identification list. Controlling FDR at the peptide-spectrum
match (PSM) level is particularly important for phosphoproteomics due to the added complexity of
modification site assignment and the lower abundance of phosphorylated species.

Link: https://openms.de/documentation/html/TOPP_FalseDiscoveryRate.html

```shell
for base_name in "${FILES[@]}"; do
    input_file="${base_name}_indexed.idXML"
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

Link: https://openms.de/documentation/html/TOPP_IDFilter.html

```shell
for base_name in "${FILES[@]}"; do
    fdr_file="${base_name}_fdr.idXML"
    filtered_file="${base_name}_filtered.idXML"

    echo "Filtering identifications..."

    # Filter identifications by score and quality
    IDFilter -in "$fdr_file" \
        -out "$filtered_file" \
        -score:peptide 0.01 \
        -remove_decoys \
        -precursor:length '8:50' \
        -precursor:charge '2:5'

    echo "IDFilter completed for: $output_file"
done
```

More stringent filtering for phospho analysis:

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

## Step 6: Phosphosite Localisation <a name="phosphositeslocalisation"></a>

Phosphosite localisation determines the exact amino acid residue(s) within a peptide that carry the
phosphorylation modification, which is critical since multiple serine, threonine, or tyrosine
residues may be present in the same peptide. Algorithms like PhosphoRS or Ascore analyse
fragmentation patterns to calculate localisation probabilities for each potential modification site.
High-confidence site localisation (typically >75% probability) is essential for understanding the
biological significance of specific phosphorylation events and their regulatory roles.

Link: https://openms.de/documentation/TOPP_PhosphoScoring.html

```shell
for base_name in "${FILES[@]}"; do
    input_file="${base_name}_picked.mzML"
    index_file="${base_name}_filtered.idXML"
    output_file="${base_name}_phosphoscoring_results.idXML"
    echo "Processing $input_file..."

    # Use PhosphoScoring for site localisation
    PhosphoScoring -in "$input_file" \
        -id "$index_file" \
        -out "$output_file" \
        -threads 4

    echo "PhosphoScoring completed for: $output_file"
done
```

## Step 7: Results Export and Analysis <a name="resultsexport"></a>

### Step 7.1: Export to Text Format and extract phosposite positions

We are simply converting the idXML file to tabular format file (mzTab), and then we use the
following python script
([ms_extract_phosphosites_from_mztab.py](../../scripts/ms_extract_phosphosites_from_mztab.py)) to
extract phosphosites from the mzTab files to a tabular file (.tsv), easier to view in spreadsheet
editor like LibreOffice.

Link: https://www.openms.org/documentation/html/TOPP_TextExporter.html

```shell
for base_name in "${FILES[@]}"; do
    input_file="${base_name}_phosphoscoring_results.idXML"
    mztab_file="${base_name}_phosphosite_results.mzTab"
    phospho_file="${base_name}_phosphosite_results.tsv"
    echo "Processing TextExporter on $input_file..."

    # Export identification results
    MzTabExporter -in "$input_file" \
        -out "$mztab_file" \
        -threads 4

    # Extract phosposite positions on protein
    ms_extract_phosphosites_from_mztab.py -in "$mztab_file" \
        -out "$phospho_file"

    echo "MzTabExporter completed for: $output_file"
done
```

## Troubleshooting

### Common Issues

1. **Low identification rate**: Check search parameters, enzyme specificity
2. **Poor localisation**: Verify fragment mass tolerance, neutral losses
3. **High FDR**: Increase database size, check decoy generation

### Performance Tips

- Use SSD storage for faster I/O
- Allocate sufficient RAM (8GB+ recommended)
- Consider parallel processing for multiple files

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
