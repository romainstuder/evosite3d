# Life-Science MCP Connectors

A reference list of MCP (Model Context Protocol) connectors for biology, pharma, biotech and cheminformatics work. It covers what each connector does, its endpoint URL, its main tools and example prompts.

> Compiled 22 October 2026 from the Claude connectors directory. Tool lists and availability may change, so check the directory for the latest details.

---

## Contents

1. [How to use these connectors](#how-to-use-these-connectors)
2. [Quick reference](#quick-reference)
3. [Cheminformatics and structural biology](#cheminformatics-and-structural-biology)
4. [Targets, genomics and bioinformatics](#targets-genomics-and-bioinformatics)
5. [Literature](#literature)
6. [Pharma and clinical](#pharma-and-clinical)
7. [Gaps and workarounds](#gaps-and-workarounds)
8. [Example workflows](#example-workflows)

---

## How to use these connectors

**In Claude (web, desktop or mobile):** open the connectors directory, find the connector by name and connect it. Connectors marked **Open** need no account. The others will ask you to sign in to the provider. Once a connector is connected and enabled for a chat, ask for what you need in plain language and Claude will call its tools.

**In other MCP clients:** each connector listed here is a remote MCP server. If your client supports remote (HTTP) MCP servers, you can add the endpoint URL as a server. Connectors that need sign-in will start an OAuth flow, and some require an account or licence with the provider.

**Tips**

- Name the database in your prompt (for example, "using ChEMBL…") when more than one connector could answer.
- Ask for identifiers in the output, such as ChEMBL IDs, PMIDs, DOIs or NCT numbers, so you can verify results.
- For large pulls, ask Claude to save the results to a file and analyse them in code rather than in the chat.

---

## Quick reference

| Connector         | Area                        | Endpoint URL                                      | Access  |
| ----------------- | --------------------------- | ------------------------------------------------- | ------- |
| ChEMBL            | Cheminformatics             | `https://hcls.mcp.claude.com/chembl/mcp`          | Open    |
| Boltz API         | Structure / binding         | `https://mcp.boltz.bio/mcp`                       | Sign-in |
| Inductive Bio     | ADMET prediction            | `https://api.inductive.bio/mcp`                   | Sign-in |
| Open Targets      | Target discovery            | `https://mcp.platform.opentargets.org/mcp`        | Open    |
| LatchBio          | Bioinformatics workflows    | `https://mcp.latch.bio/mcp`                       | Sign-in |
| Synthesize Bio    | Gene expression             | `https://app.synthesize.bio/api/mcp`              | Sign-in |
| Helix GenoSphere  | Human genomics              | `https://mcp.hrn-production.helix.com/mcp`        | Sign-in |
| Owkin             | AI agents for biology       | `https://mcp.k.owkin.com/mcp`                     | Sign-in |
| PubMed            | Literature                  | `https://pubmed.mcp.claude.com/mcp`               | Open    |
| bioRxiv / medRxiv | Preprints                   | `https://hcls.mcp.claude.com/biorxiv/mcp`         | Open    |
| Scite             | Citation evidence           | `https://api.scite.ai/mcp`                        | Sign-in |
| Clinical Trials   | ClinicalTrials.gov          | `https://hcls.mcp.claude.com/clinical_trials/mcp` | Open    |
| AdisInsight       | Drug and trial intelligence | `https://adisinsight-mcp.springer.com/mcp`        | Sign-in |
| Amass             | Unified life-science data   | `https://mcp.amass.tech/mcp`                      | Sign-in |
| Medidata          | Trial site ranking          | Directory only                                    | Sign-in |
| Benchling         | ELN / registry              | Directory only                                    | Sign-in |

"Directory only" means no public endpoint is listed, so connect through the Claude connectors directory.

---

## Cheminformatics and structural biology

### ChEMBL

Access to EMBL-EBI's database of bioactive molecules.

- **URL:** `https://hcls.mcp.claude.com/chembl/mcp`
- **Access:** Open
- **Tools:** `compound_search`, `drug_search`, `target_search`, `get_bioactivity`, `get_mechanism`, `get_admet`
- **Example prompts:**
  - "Find ChEMBL compounds with IC50 < 100 nM against EGFR and list their SMILES."
  - "What is the mechanism of action of imatinib according to ChEMBL?"

### Boltz API

Biomolecular structure and binding prediction, virtual screening and binder design.

- **URL:** `https://mcp.boltz.bio/mcp`
- **Access:** Sign-in (compute jobs may be billed)
- **Tools:** `boltz_estimate_structure_and_binding`, `boltz_get_structure_and_binding_prediction`, `boltz_estimate_protein_screen`, `boltz_start_protein_screen`, `boltz_estimate_small_molecule_screen`, `boltz_get_screen_results`, `boltz_get_job_status`, `boltz_get_guidance`
- **Example prompts:**
  - "Estimate the cost of predicting the structure of this sequence bound to this ligand SMILES, then run it."
  - "Screen this 500-compound library against my target and rank the hits by predicted affinity."
- **Note:** jobs run asynchronously. Use the estimate tools first, then poll job status.

### Inductive Bio

ADMET property prediction models for drug discovery.

- **URL:** `https://api.inductive.bio/mcp`
- **Access:** Sign-in
- **Tools:** `list_available_models`, `predict_properties`
- **Example prompts:**
  - "Which ADMET models are available?"
  - "Predict solubility, permeability and hERG liability for these SMILES."

---

## Targets, genomics and bioinformatics

### Open Targets

Evidence linking targets to diseases, used for target identification and prioritisation.

- **URL:** `https://mcp.platform.opentargets.org/mcp`
- **Access:** Open
- **Tools:** `search_entities`, `get_open_targets_graphql_schema`, `query_open_targets_graphql`, `batch_query_open_targets_graphql`
- **Example prompts:**
  - "Top 20 targets associated with idiopathic pulmonary fibrosis, with their overall association scores."
  - "Which approved drugs target genes associated with type 2 diabetes?"

### LatchBio

Run bioinformatics workflows and manage data on the Latch platform.

- **URL:** `https://mcp.latch.bio/mcp`
- **Access:** Sign-in (Latch account)
- **Tools:** `list_files`, `get_file`, `list_workflows`, `get_workflow_schema`, `launch_workflow`, `list_executions`, `get_execution`, `get_task_logs`
- **Example prompts:**
  - "List my workflows and launch the RNA-seq pipeline on the FASTQs in /project/run42."
  - "Why did my last execution fail? Show the task logs."

### Synthesize Bio

Generates and analyses gene expression from a "virtual human" model.

- **URL:** `https://app.synthesize.bio/api/mcp`
- **Access:** Sign-in
- **Tools:** `resolve_sample_metadata`, `analyze_gene_expression`, `get_analysis_results`, `get_saved_analysis_report`, `annotate_genes`, `get_counts_data_url`
- **Example prompt:** "Simulate expression in lung tissue with and without the KRAS G12C condition and run a differential expression analysis."

### Helix GenoSphere

Human genomic data linked to longitudinal clinical records.

- **URL:** `https://mcp.hrn-production.helix.com/mcp`
- **Access:** Sign-in (Helix access required)
- **Tools:** `get_dataset_overview`, `check_gene_availability`, `get_variant_frequency`, `search_concepts`, `get_patient_count`
- **Example prompt:** "What is the frequency of this PCSK9 variant, and how many carriers have a hypercholesterolaemia diagnosis?"

### Owkin

AI agents for biology, focused on histopathology and cohort analysis.

- **URL:** `https://mcp.k.owkin.com/mcp`
- **Access:** Sign-in
- **Tools:** `list_cohorts`, `list_cell_types`, `filter_slides`, `get_histomics`, `display_slide_thumbnail`, `survival_analysis`, `features_description`, `histoplus_help`
- **Example prompt:** "Run a survival analysis on the cohort, stratified by tumour-infiltrating lymphocyte density."

---

## Literature

### PubMed

Search and retrieve biomedical literature.

- **URL:** `https://pubmed.mcp.claude.com/mcp`
- **Access:** Open
- **Tools:** `search_articles`, `get_article_metadata`, `find_related_articles`, `lookup_article_by_citation`, `convert_article_ids`, `get_full_text_article`, `get_copyright_status`
- **Example prompts:**
  - "Find papers since 2020 on RubisCO stability–activity trade-offs."
  - "Get articles related to PMID 24469821."

### bioRxiv / medRxiv

Preprint search and metadata.

- **URL:** `https://hcls.mcp.claude.com/biorxiv/mcp`
- **Access:** Open
- **Tools:** `search_preprints`, `search_biorxiv_publications`, `get_preprint`, `get_categories`, `search_by_funder`, `search_published_articles`, `search_publisher_articles`, `get_content_statistics`
- **Example prompt:** "Recent bioRxiv preprints on ancestral sequence reconstruction of enzymes."

### Scite

Literature search that shows whether citing papers support or contradict a claim.

- **URL:** `https://api.scite.ai/mcp`
- **Access:** Sign-in
- **Tools:** `search_literature`
- **Example prompt:** "Is there evidence supporting or contradicting the claim that M309I is a C4 catalytic switch in RubisCO?"

---

## Pharma and clinical

### Clinical Trials

Data from ClinicalTrials.gov.

- **URL:** `https://hcls.mcp.claude.com/clinical_trials/mcp`
- **Access:** Open
- **Tools:** `search_trials`, `get_trial_details`, `search_by_eligibility`, `search_by_sponsor`, `search_investigators`, `analyze_endpoints`
- **Example prompt:** "Phase 3 trials currently recruiting for NASH. Summarise their primary endpoints."

### AdisInsight (Springer Nature)

Drug development pipeline and clinical trial intelligence.

- **URL:** `https://adisinsight-mcp.springer.com/mcp`
- **Access:** Sign-in (subscription)
- **Tools:** `search_drugs`, `get_drug`, `search_trials`, `get_trial`, `search_drug_companies`, `search_trial_companies`, `generate_chart`
- **Example prompt:** "Pipeline of GLP-1 agonists by development phase, with a chart."

### Amass

A single cross-linked source covering publications, trials, drugs, genes, regulatory approvals and patents.

- **URL:** `https://mcp.amass.tech/mcp`
- **Access:** Sign-in
- **Tools:** `search_amass_biomedcore_records`, `get_amass_biomedcore_record`, `get_amass_drugcore_record`, `get_amass_genecore_record`, `get_amass_trialcore_record`, `get_amass_patentcore_record`, `get_amass_regulatorycore_record`, `get_amass_regulatorycore_document_section`, and others
- **Example prompt:** "For tirzepatide, link the key trials, regulatory approvals and core patents."

### Medidata

Clinical trial software and site-ranking tools.

- **URL:** Directory only
- **Access:** Sign-in (Medidata account)
- **Tools:** `get_ranked_sites`, `get_supported_disease_areas`, `search_support_documentation`
- **Example prompt:** "Rank sites for a phase 2 oncology trial in non-small cell lung cancer."

### Benchling

Electronic lab notebook and registry, with Benchling's scientific AI.

- **URL:** Directory only
- **Access:** Sign-in (Benchling tenant)
- **Tools:** `benchling_agent`, `benchling_task_status`, `benchling_agent_edit_entry`
- **Example prompt:** "Summarise last week's notebook entries for the RubisCO expression project."

---

## Gaps and workarounds

There is currently no dedicated connector for **UniProt**, **RCSB PDB**, **Ensembl** or **PubChem**. You have two options:

- Ask Claude to fetch the public files or API results directly (for example, PDB `1WDD`, rice RubisCO) and analyse them in Python with Biopython, RDKit or similar libraries.
- Use **Amass** or **Open Targets** for gene and protein annotation, and **ChEMBL** for compound data.

---

## Example workflows

**Structure-based enzyme evolution study** (for example, RubisCO C3→C4 stability analysis)

1. **PubMed / bioRxiv:** gather the literature on RubisCO trade-offs and ancestral reconstruction.
2. **Boltz:** predict structures for ancestral or mutant sequences.
3. Analyse in code: compute residue contacts and interfaces, and compare models.
4. **Scite:** check how later work supports or challenges the key claims.

**Early drug discovery triage**

1. **Open Targets:** shortlist targets for a disease.
2. **ChEMBL:** pull known actives and their bioactivity for each target.
3. **Inductive Bio:** predict ADMET properties for the candidate compounds.
4. **Clinical Trials / AdisInsight:** check the competitive landscape.
