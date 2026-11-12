# Life-Science MCP Connectors

A reference list of Model Context Protocol (MCP) connectors for biology, pharma, biotech and cheminformatics work. It covers what each connector does, its endpoint URL, its main tools and example prompts.

> Compiled 22 October 2026 from the Claude connectors directory and provider announcements. Tool lists and availability may change, so check the directory for the latest details.

---

## Contents

1. [How to use these connectors](#how-to-use-these-connectors)
2. [Quick reference](#quick-reference)
3. [Cheminformatics and structural biology](#cheminformatics-and-structural-biology)
4. [Targets, genomics and bioinformatics](#targets-genomics-and-bioinformatics)
5. [Literature](#literature)
6. [Pharma and clinical](#pharma-and-clinical)
7. [Figures and illustration](#figures-and-illustration)
8. [Other directories](#other-directories)
9. [Gaps and workarounds](#gaps-and-workarounds)
10. [Example workflows](#example-workflows)

---

## How to use these connectors

**In Claude (web, desktop or mobile):** open the connectors directory, find the connector by name and connect it. Connectors marked **Open** need no account. The others will ask you to sign in to the provider, and some need a paid plan or a customer relationship (see the access labels under Quick reference). Once a connector is connected and enabled for a chat, ask for what you need in plain language and Claude will call its tools.

**In other MCP clients:** each connector listed here is a remote MCP server. If your client supports remote (HTTP) MCP servers, you can add the endpoint URL as a server. Connectors that need sign-in will start an OAuth flow, and some require an account or licence with the provider.

**Tips**

- Name the database in your prompt (for example, "using ChEMBL…") when more than one connector could answer.
- Ask for identifiers in the output, such as ChEMBL IDs, PMIDs, DOIs or NCT numbers, so you can verify results.
- For large pulls, ask Claude to save the results to a file and analyse them in code rather than in the chat.

---

## Quick reference

### Open connectors

No account needed. Add each one by its endpoint URL, or find it in the Claude connectors directory.

| Connector         | Area               | Endpoint URL                                      |
| ----------------- | ------------------ | ------------------------------------------------- |
| ChEMBL            | Cheminformatics    | `https://hcls.mcp.claude.com/chembl/mcp`          |
| RCSB PDB (beta)   | Structure archive  | `https://mcp-beta.rcsb.org/mcp`                   |
| Open Targets      | Target discovery   | `https://mcp.platform.opentargets.org/mcp`        |
| PubMed            | Literature         | `https://pubmed.mcp.claude.com/mcp`               |
| bioRxiv / medRxiv | Preprints          | `https://hcls.mcp.claude.com/biorxiv/mcp`         |
| Clinical Trials   | ClinicalTrials.gov | `https://hcls.mcp.claude.com/clinical_trials/mcp` |

### Built into Claude Science

Featured connectors that are on by default in Claude Science. They are read-only, need no account or key, and have no public endpoint listed. Claude Science also includes many more database connectors, such as Ensembl, UniProt and PubChem. See the [Claude Science connectors page](https://claude.com/docs/claude-science/connectors-and-skills) for the full list.

| Connector | Area                  |
| --------- | --------------------- |
| BioMart   | Gene annotation       |
| CellGuide | Cell types            |
| ZINC      | Purchasable chemistry |

### Connectors that need an account

| Connector        | Area                             | Endpoint URL                               | Access                   |
| ---------------- | -------------------------------- | ------------------------------------------ | ------------------------ |
| Boltz API        | Structure / binding              | `https://mcp.boltz.bio/mcp`                | Pay-per-use              |
| Inductive Bio    | ADMET prediction                 | `https://api.inductive.bio/mcp`            | Account, pricing unclear |
| LatchBio         | Bioinformatics workflows         | `https://mcp.latch.bio/mcp`                | Pay-per-use              |
| Synthesize Bio   | Gene expression                  | `https://app.synthesize.bio/api/mcp`       | Account, pricing unclear |
| Helix GenoSphere | Human genomics                   | `https://mcp.hrn-production.helix.com/mcp` | Partners only            |
| Owkin            | AI agents for biology            | `https://mcp.k.owkin.com/mcp`              | Account, pricing unclear |
| Scite            | Citation evidence                | `https://api.scite.ai/mcp`                 | Account, pricing unclear |
| AdisInsight      | Drug and trial intelligence      | `https://adisinsight-mcp.springer.com/mcp` | Subscription             |
| Amass            | Unified life-science data        | `https://mcp.amass.tech/mcp`               | Account, pricing unclear |
| Medidata         | Trial site ranking               | Directory only                             | Customers only           |
| Benchling        | ELN / registry                   | Directory only                             | Customers only           |
| BioRender        | Scientific figures               | `https://mcp.services.biorender.com/mcp`   | Pay-per-use              |
| 10x Genomics     | Single-cell and spatial analysis | Desktop extension                          | Account, pricing unclear |
| Synapse.org      | Research data repository         | `https://mcp.synapse.org/mcp`              | Account, pricing unclear |
| Scholar Gateway  | Wiley literature                 | `https://connector.scholargateway.ai/mcp`  | Account, pricing unclear |
| Consensus        | Literature search                | `https://mcp.consensus.app/mcp`            | Free tier, paid extras   |

"Directory only" means no public endpoint is listed, so connect through the Claude connectors directory. "Desktop extension" means the connector is installed as a Claude Desktop extension, not added by URL.

**Access labels** (for connectors that need an account)

- **Pay-per-use:** needs an account, and usage is billed.
- **Subscription:** needs a paid subscription with the provider.
- **Customers only / Partners only:** limited to the provider's customers or approved partners.
- **Free tier, paid extras:** usable without paying, with higher limits on paid plans.
- **Account, pricing unclear:** needs an account with the provider. No pricing or eligibility information was found, so check with the provider.

These labels come from the providers' pages and third-party listings, not from the providers' pricing pages. Confirm costs before relying on a connector.

---

## Cheminformatics and structural biology

### ChEMBL

Access to EMBL-EBI's database of bioactive molecules.

- **Website URL:** <https://www.ebi.ac.uk/chembl/>
- **Documentation URL:** <https://academy.claude.com/tutorials/using-the-chembl-connector-in-claude>
- **Endpoint URL:** `https://hcls.mcp.claude.com/chembl/mcp`
- **Access:** Open. No account needed.
- **Tools:** `compound_search`, `drug_search`, `target_search`, `get_bioactivity`, `get_mechanism`, `get_admet`
- **Example prompts:**
  - "Find ChEMBL compounds with IC50 < 100 nM against EGFR and list their SMILES."
  - "What is the mechanism of action of imatinib according to ChEMBL?"

### RCSB PDB (beta)

Search and retrieve experimentally determined 3D structures from the Protein Data Bank. Open source: [rcsb/rcsb-mcp](https://github.com/rcsb/rcsb-mcp).

- **Website URL:** <https://www.rcsb.org>
- **Documentation URL:** <https://github.com/rcsb/rcsb-mcp>
- **Endpoint URL:** `https://mcp-beta.rcsb.org/mcp`
- **Access:** Open. The announcement mentions no sign-in.
- **Capabilities:** the announcement groups them by API rather than listing tool names.
  - Search API: find structures by keyword, attributes, sequence similarity, chemical structure, 3D shape, or sequence and structural motifs, combined with Boolean logic and summarised with counts and facets.
  - Data API: details on entries, entities, assemblies, chains, ligands and interfaces, with linked UniProt and PubMed records.
  - Sequence Coordinates API: map sequences and positional features between PDB, UniProt and NCBI.
  - Helper tools: convert free text into identifiers for Gene Ontology, InterPro/Pfam, EC numbers, MONDO disease terms and NCBI Taxonomy.
- **Example prompts:**
  - "Is there any protein interacting with uranium compounds?"
  - "List the residue binding sites for those results and the most frequent amino acid types."
- **Note:** this is a beta release. Tools and behaviour may change, so check key results on RCSB.org. Feedback goes to the [issue tracker](https://github.com/rcsb/rcsb-mcp/issues). Source: [RCSB announcement](https://www.rcsb.org/news/6acfa3d74ec37840b7ddf62e).

### Boltz API

Biomolecular structure and binding prediction, virtual screening and binder design.

- **Website URL:** <https://boltz.bio>
- **Documentation URL:** <https://claude.com/marketplace/connectors/boltz-api>
- **Endpoint URL:** `https://mcp.boltz.bio/mcp`
- **Access:** Pay-per-use. Sign in to a Boltz account; compute jobs may be billed. Use the estimate tools first.
- **Tools:** `boltz_estimate_structure_and_binding`, `boltz_get_structure_and_binding_prediction`, `boltz_estimate_protein_screen`, `boltz_start_protein_screen`, `boltz_estimate_small_molecule_screen`, `boltz_get_screen_results`, `boltz_get_job_status`, `boltz_get_guidance`
- **Example prompts:**
  - "Estimate the cost of predicting the structure of this sequence bound to this ligand SMILES, then run it."
  - "Screen this 500-compound library against my target and rank the hits by predicted affinity."
- **Note:** jobs run asynchronously. Use the estimate tools first, then poll job status.

### Inductive Bio

ADMET property prediction models for drug discovery.

- **Website URL:** <https://www.inductive.bio>
- **Documentation URL:** none. Inductive Bio has no documentation page for its MCP server. The
  [Claude connector listing](https://claude.com/marketplace/connectors/inductive-bio) is the closest reference.
- **Endpoint URL:** `https://api.inductive.bio/mcp`
- **Access:** Account, pricing unclear. Sign in to an Inductive Bio account; no pricing or eligibility information found.
- **Tools:** `list_available_models`, `predict_properties`
- **Example prompts:**
  - "Which ADMET models are available?"
  - "Predict solubility, permeability and hERG liability for these SMILES."

### ZINC

Purchasable chemical space, for finding compounds that can be bought. A Featured connector in Claude Science.

- **Website URL:** <https://zinc.docking.org>
- **Documentation URL:** <https://claude.com/docs/claude-science/connectors-and-skills>
- **Endpoint URL:** Built into Claude Science (no public endpoint listed)
- **Access:** Open. Featured connectors are read-only and need no account or key. Some underlying databases have non-commercial or attribution terms, so check each source's licence.
- **Tools:** not listed in the documentation.

---

## Targets, genomics and bioinformatics

### Open Targets

Evidence linking targets to diseases, used for target identification and prioritisation.

- **Website URL:** <https://platform.opentargets.org>
- **Documentation URL:** <https://platform-docs.opentargets.org/data-access/model-context-protocol>
- **Endpoint URL:** `https://mcp.platform.opentargets.org/mcp`
- **Access:** Open. No account needed.
- **Tools:** `search_entities`, `get_open_targets_graphql_schema`, `query_open_targets_graphql`, `batch_query_open_targets_graphql`
- **Example prompts:**
  - "Top 20 targets associated with idiopathic pulmonary fibrosis, with their overall association scores."
  - "Which approved drugs target genes associated with type 2 diabetes?"

### LatchBio

Run bioinformatics workflows and manage data on the Latch platform.

- **Website URL:** <https://latch.bio>
- **Documentation URL:** <https://wiki.latch.bio/agent/latch-mcp>
- **Endpoint URL:** `https://mcp.latch.bio/mcp`
- **Access:** Pay-per-use. Sign in to a Latch account; actions run through MCP are charged as usual, and launching a workflow bills the workspace.
- **Tools:** `list_files`, `get_file`, `list_workflows`, `get_workflow_schema`, `launch_workflow`, `list_executions`, `get_execution`, `get_task_logs`
- **Example prompts:**
  - "List my workflows and launch the RNA-seq pipeline on the FASTQs in /project/run42."
  - "Why did my last execution fail? Show the task logs."

### Synthesize Bio

Generates and analyses gene expression from a "virtual human" model.

- **Website URL:** <https://www.synthesize.bio>
- **Documentation URL:** <https://docs.synthesize.bio/platform>
- **Endpoint URL:** `https://app.synthesize.bio/api/mcp`
- **Access:** Account, pricing unclear. Sign in to a Synthesize Bio account; no pricing or eligibility information found.
- **Tools:** `resolve_sample_metadata`, `analyze_gene_expression`, `get_analysis_results`, `get_saved_analysis_report`, `annotate_genes`, `get_counts_data_url`
- **Example prompt:** "Simulate expression in lung tissue with and without the KRAS G12C condition and run a differential expression analysis."

### Helix GenoSphere

Human genomic data linked to longitudinal clinical records.

- **Website URL:** <https://www.helix.com>
- **Documentation URL:** none. Helix has no documentation page for its MCP server. The
  [launch press release](https://www.helix.com/press-releases/helix-genosphere-mcp-connector-human-genetics) is the
  closest reference.
- **Endpoint URL:** `https://mcp.hrn-production.helix.com/mcp`
- **Access:** Partners only. Limited to approved enterprise health-system and life-sciences partners of Helix.
- **Tools:** `get_dataset_overview`, `check_gene_availability`, `get_variant_frequency`, `search_concepts`, `get_patient_count`
- **Example prompt:** "What is the frequency of this PCSK9 variant, and how many carriers have a hypercholesterolaemia diagnosis?"

### Owkin

AI agents for biology, focused on histopathology and cohort analysis.

- **Website URL:** <https://www.owkin.com>
- **Documentation URL:** <https://docs.owkin.com/connect-and-integrate/mcp-and-agent-tools>
- **Endpoint URL:** `https://mcp.k.owkin.com/mcp`
- **Access:** Account, pricing unclear. Sign in to an Owkin account; no pricing or eligibility information found.
- **Tools:** `list_cohorts`, `list_cell_types`, `filter_slides`, `get_histomics`, `display_slide_thumbnail`, `survival_analysis`, `features_description`, `histoplus_help`
- **Example prompt:** "Run a survival analysis on the cohort, stratified by tumour-infiltrating lymphocyte density."

### 10x Genomics

Run single-cell and spatial analyses on 10x Cloud Analysis from a prompt.

- **Website URL:** <https://www.10xgenomics.com>
- **Documentation URL:** <https://support.claude.com/en/articles/12614803-using-the-10x-genomics-extension-in-claude>
- **Endpoint URL:** none. This is a Claude Desktop extension, not a remote URL.
- **Access:** Account, pricing unclear. Needs a 10x Cloud Analysis account and a 10x Cloud access token (Account Settings, Security).
- **Capabilities:** set up runs, launch pipelines, monitor status and download results, covering alignment, Feature Barcode matrix generation and clustering. Exact tool names were not found.
- **Example prompt:** "Align these FASTQs, generate the feature barcode matrix and cluster the cells."

### Synapse.org

Sage Bionetworks' research data repository, for finding and reading datasets and their annotations.

- **Website URL:** <https://www.synapse.org>
- **Documentation URL:** <https://github.com/Sage-Bionetworks/synapse-mcp>
- **Endpoint URL:** `https://mcp.synapse.org/mcp`
- **Access:** Account, pricing unclear. Sign in with a Synapse account (OAuth).
- **Tools:** `get_entity`, `get_entity_annotations`, `get_entity_provenance`, `get_entity_children`, `search_synapse`
- **Example prompt:** "Search Synapse for Alzheimer's RNA-seq datasets and show the provenance of the first result."
- **Note:** the directory lists these five read-only tools. The repository describes a different default, with `get_entity` and `search_synapse` always visible and other tools reached through search, so the toolset may vary by client or version. Not to be confused with SynapseNetwork, an unrelated product.

### BioMart

Query BioMart marts for gene, transcript and variant annotation. A Featured connector in Claude Science.

- **Website URL:** <https://www.biomart.org>
- **Documentation URL:** <https://claude.com/docs/claude-science/connectors-and-skills>
- **Endpoint URL:** Built into Claude Science (no public endpoint listed)
- **Access:** Open. Featured connectors are read-only and need no account or key.
- **Tools:** not listed in the documentation.
- **Note:** a community server, [jzinno/biomart-mcp](https://github.com/jzinno/biomart-mcp), runs locally for other clients.

### CellGuide

Cell-type information from CZ CELLxGENE CellGuide. A Featured connector in Claude Science.

- **Website URL:** <https://cellxgene.cziscience.com/cellguide>
- **Documentation URL:** <https://claude.com/docs/claude-science/connectors-and-skills>
- **Endpoint URL:** Built into Claude Science (no public endpoint listed)
- **Access:** Open. Featured connectors are read-only and need no account or key.
- **Tools:** not listed in the documentation.

---

## Literature

### PubMed

Search and retrieve biomedical literature.

- **Website URL:** <https://pubmed.ncbi.nlm.nih.gov>
- **Documentation URL:** <https://academy.claude.com/tutorials/using-the-pubmed-connector-in-claude>
- **Endpoint URL:** `https://pubmed.mcp.claude.com/mcp`
- **Access:** Open. No account needed.
- **Tools:** `search_articles`, `get_article_metadata`, `find_related_articles`, `lookup_article_by_citation`, `convert_article_ids`, `get_full_text_article`, `get_copyright_status`
- **Example prompts:**
  - "Find papers since 2020 on RubisCO stability–activity trade-offs."
  - "Get articles related to PMID 24469821."

### bioRxiv / medRxiv

Preprint search and metadata.

- **Website URL:** <https://www.biorxiv.org>
- **Documentation URL:** <https://academy.claude.com/tutorials/using-the-biorxiv-and-medrxiv-connector-in-claude>
- **Endpoint URL:** `https://hcls.mcp.claude.com/biorxiv/mcp`
- **Access:** Open. No account needed.
- **Tools:** `search_preprints`, `search_biorxiv_publications`, `get_preprint`, `get_categories`, `search_by_funder`, `search_published_articles`, `search_publisher_articles`, `get_content_statistics`
- **Example prompt:** "Recent bioRxiv preprints on ancestral sequence reconstruction of enzymes."

### Scite

Literature search that shows whether citing papers support or contradict a claim.

- **Website URL:** <https://scite.ai>
- **Documentation URL:** <https://docs.scite.ai/mcp/overview>
- **Endpoint URL:** `https://api.scite.ai/mcp`
- **Access:** Account, pricing unclear. Sign in to a Scite account; sources disagree on whether a paid plan or free trial is required.
- **Tools:** `search_literature`
- **Example prompt:** "Is there evidence supporting or contradicting the claim that M309I is a C4 catalytic switch in RubisCO?"

### Scholar Gateway

Semantic search over Wiley's scholarly content, returning passages with citation metadata.

- **Website URL:** <https://www.wiley.com/en-us/research/products/wiley-scholar-gateway>
- **Documentation URL:** <https://academy.claude.com/tutorials/using-the-scholar-gateway-connector-in-claude>
- **Endpoint URL:** `https://connector.scholargateway.ai/mcp`
- **Access:** Account, pricing unclear. Sign in with a Wiley account. Wiley's page advertises a free tier of 30 queries a month, while the Claude Academy guide says existing Wiley journal subscribers will need to upgrade after a trial period.
- **Tools:** a single semantic search tool (community-reported name: `semanticSearch`)
- **Example prompt:** "Find Wiley papers on thermostable RubisCO variants and quote the relevant passages."
- **Note:** covers Wiley content only.

### Consensus

Search for scientific papers, with answers grounded in the literature.

- **Website URL:** <https://consensus.app>
- **Documentation URL:** <https://docs.consensus.app/docs/mcp>
- **Endpoint URL:** `https://mcp.consensus.app/mcp`
- **Access:** Free tier, paid extras. It works without signing in at reduced limits. Paid plans can enable extra usage at about $0.05 per call, according to the docs.
- **Tools:** search. The docs mention five tools, but a fetch tool is currently available for ChatGPT only, so Claude users may see only search.
- **Example prompt:** "What does the evidence say about ancestral sequence reconstruction improving enzyme thermostability?"
- **Note:** the URL must end in `/mcp`. A call covers up to 100 papers, and the rate limit is one request per second.

---

## Pharma and clinical

### Clinical Trials

Data from ClinicalTrials.gov.

- **Website URL:** <https://clinicaltrials.gov>
- **Documentation URL:** <https://academy.claude.com/tutorials/using-the-clinicaltrials-gov-connector-in-claude>
- **Endpoint URL:** `https://hcls.mcp.claude.com/clinical_trials/mcp`
- **Access:** Open. No account needed.
- **Tools:** `search_trials`, `get_trial_details`, `search_by_eligibility`, `search_by_sponsor`, `search_investigators`, `analyze_endpoints`
- **Example prompt:** "Phase 3 trials currently recruiting for MASH (formerly NASH). Search both terms and summarise their primary endpoints."

### AdisInsight (Springer Nature)

Drug development pipeline and clinical trial intelligence.

- **Website URL:** <https://adisinsight.springer.com>
- **Documentation URL:** <https://adisinsight-mcp.springer.com/>
- **Endpoint URL:** `https://adisinsight-mcp.springer.com/mcp`
- **Access:** Subscription. Requires an AdisInsight subscription and a sign-in.
- **Tools:** `search_drugs`, `get_drug`, `search_trials`, `get_trial`, `search_drug_companies`, `search_trial_companies`, `generate_chart`
- **Example prompt:** "Pipeline of GLP-1 agonists by development phase, with a chart."

### Amass

A single cross-linked source covering publications, trials, drugs, genes, regulatory approvals and patents.

- **Website URL:** <https://amass.tech>
- **Documentation URL:** <https://amass.tech/mcp>
- **Endpoint URL:** `https://mcp.amass.tech/mcp`
- **Access:** Account, pricing unclear. Sign in to an Amass account (OAuth), or use an Amass API key as a bearer token; no pricing information found.
- **Tools:** `search_amass_biomedcore_records`, `get_amass_biomedcore_record`, `get_amass_drugcore_record`, `get_amass_genecore_record`, `get_amass_trialcore_record`, `get_amass_patentcore_record`, `get_amass_regulatorycore_record`, `get_amass_regulatorycore_document_section`, and others
- **Example prompt:** "For tirzepatide, link the key trials, regulatory approvals and core patents."

### Medidata

Clinical trial software and site-ranking tools.

- **Website URL:** <https://www.medidata.com>
- **Documentation URL:** <https://academy.claude.com/tutorials/using-the-medidata-connector-in-claude>
- **Endpoint URL:** Directory only
- **Access:** Customers only. For Medidata customers; site ranking needs Intelligent Trials.
- **Tools:** `get_ranked_sites`, `get_supported_disease_areas`, `search_support_documentation`
- **Example prompt:** "Rank sites for a phase 2 oncology trial in non-small cell lung cancer."

### Benchling

Electronic lab notebook and registry, with Benchling's scientific AI.

- **Website URL:** <https://www.benchling.com>
- **Documentation URL:** <https://help.benchling.com/hc/en-us/articles/40342713479437-Configure-Benchling-s-MCP-Server-for-other-MCP-clients>
- **Endpoint URL:** Directory only
- **Access:** Customers only. Needs a Benchling tenant with Benchling AI and the remote MCP enabled.
- **Tools:** `benchling_agent`, `benchling_task_status`, `benchling_agent_edit_entry`
- **Example prompt:** "Summarise last week's notebook entries for the RubisCO expression project."

---

## Figures and illustration

### BioRender

Search BioRender's template library and your own files, and create scientific figures from a text description.

- **Website URL:** <https://www.biorender.com>
- **Documentation URL:** <https://support.claude.com/en/articles/12614795-using-the-biorender-connector-in-claude>
- **Endpoint URL:** `https://mcp.services.biorender.com/mcp`
- **Access:** Pay-per-use. Sign in with a BioRender account (OAuth). Generating custom figure previews uses BioRender AI credits.
- **Capabilities:** search templates, search your own and shared files, and generate custom figures from a description, optionally with a reference image or template. Exact tool names were not found.
- **Example prompt:** "Create a figure of the RubisCO holoenzyme showing the large and small subunits."
- **Note:** searches, prompts and generated figures are shared with the AI assistant you connect.

---

## Other directories

This list is a short, curated set of connectors. Larger directories exist, and they also list community-built servers for databases without an official connector, such as UniProt, Ensembl and PubChem.

- [mcp.so, bioinformatics tag](https://mcp.so/tags/bioinformatics): a general MCP directory. The bioinformatics tag lists about 20 servers and clients (for example BioMCP, ENCODE, UniProt and Seqera). It does not mark which are official, and most look like community projects.
- [Glama, Biology & Medicine](https://glama.ai/mcp/categories/biology-and-medicine): Glama's category for life-science and medical servers (genomics, protein structure, clinical vocabularies, health records). Many entries carry an "official" tag, which the page does not define.
- [IntuitionLabs, Top MCP Servers for Biotech](https://intuitionlabs.ai/pdfs/top-mcp-servers-for-biotech-connecting-ai-to-research-data.pdf): a PDF report from August 2026 covering literature, genomics and proteins (NCBI Datasets, Ensembl, UniProt, AlphaFold), pathways (Reactome, STRING), chemistry and drugs (ChEMBL, PubChem, OpenFDA) and ClinicalTrials.gov.

Servers in these directories are not vetted here. Check who maintains a server, what data it sends where, and its licence before connecting it.

---

## Gaps and workarounds

There is currently no dedicated connector in the Claude connectors directory for **UniProt**, **Ensembl** or **PubChem**. In Claude Science they are built-in Featured connectors, along with GO, Reactome, ClinVar, gnomAD, AlphaFold and others, listed on the [Claude Science connectors page](https://claude.com/docs/claude-science/connectors-and-skills). RCSB PDB has a beta connector (see above), and the [other directories](#other-directories) list community servers for some of these databases. Outside Claude Science you have these options:

- Ask Claude to fetch the public files or API results directly (for example, PDB `1WDD`, rice RubisCO) and analyse them in Python with Biopython, RDKit or similar libraries.
- Use **Amass** or **Open Targets** for gene and protein annotation, and **ChEMBL** for compound data.
- Try a community server from the other directories, after checking who maintains it.

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
