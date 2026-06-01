#!/usr/bin/env nextflow

/*
 * Positive selection site-model pipeline
 *
 * Runs the CodeML site-model workflow:
 *   1a. ensembl_msa_tree_extract.py              — Bronze: raw Ensembl tree + MSA + CDS
 *   1b. fetch_structure.py                       — 3D structure (AlphaFold / RCSB)
 *   2a. ensembl_msa_tree_transform.py            — Silver: cleaned, trimmed, PHYLIP-ready
 *   2b1. run_codeml_site_models.py --model M8     — alternative model
 *   2b2. run_codeml_site_models.py --model M8a    — null model
 *   2b3. run_codeml_site_models.py --extract-only — LRT + BEB extraction
 *   2c. run_hyphy_site_models.py                 — FUBAR (parallel per-site cross-check)
 *   2d. run_codeml_branch_site_models.py         — branch-site model A along the target lineage
 *   3.  analyse_site_models.py                   — Jalview, PyMOL, dN/dS plot
 *
 * Steps 1a (Bronze) and 1b (structure) run in parallel; 2a (Silver) consumes
 * 1a; 2b1/2b2/2c run in parallel; 2b3 joins 2b1+2b2 then feeds 3.
 *
 * Usage:
 *   nextflow run main.nf --gene_symbols "HLA-DQB1,TRIM5,FOXP2"
 *   nextflow run main.nf --gene_symbols "HLA-DQB1" --uniprot P01920 --pdb 1UVQ --chain B --resi_offset 32
 *   nextflow run main.nf --input genes.csv
 */

nextflow.enable.dsl=2

// ---------------------------------------------------------------------------
//  Parameters
// ---------------------------------------------------------------------------

params.gene_symbols    = null       // Comma-separated gene symbols
params.input           = null       // CSV: gene_symbol[,uniprot,pdb,chain,resi_offset]
params.taxon           = 1437010       // NCBI taxon ID (default Boreoeutheria)
params.plot_threshold  = 0.50       // BEB probability threshold for dN/dS plot
params.outdir          = "results"
// Top-level output directory
params.codeml_method   = "M8"       // CODEML MODEL (set "" to skip)
params.skip_codeml_run = false      // Skip 2b1/2b2; reuse existing .mlc in silver
params.hyphy_method    = ""    // FUBAR, BUSTED, or BOTH (set "" to skip)
params.branch_site     = false      // Run branch-site model A along the target lineage
params.max_pdb         = 5          // Max RCSB PDBs from SIFTS (0 to disable)

// Per-gene overrides (single gene via --gene_symbols only)
params.uniprot         = null
params.pdb             = null
params.chain           = "A"
params.resi_offset     = 0

// Path to Python scripts (auto-detected)
params.script_dir      = "${projectDir}"


// ---------------------------------------------------------------------------
//  Input channel
// ---------------------------------------------------------------------------

def build_input_channel() {
    if (params.input) {
        return Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row ->
                tuple(
                    row.gene_symbol,
                    row.uniprot    ?: '',
                    row.pdb        ?: '',
                    row.chain      ?: 'A',
                    (row.resi_offset ?: '0') as int
                )
            }
    } else if (params.gene_symbols) {
        return Channel
            .of(params.gene_symbols.tokenize(','))
            .flatten()
            .map { symbol ->
                tuple(
                    symbol.trim(),
                    params.uniprot      ?: '',
                    params.pdb          ?: '',
                    params.chain        ?: 'A',
                    (params.resi_offset ?: 0) as int
                )
            }
    } else {
        error "Provide --gene_symbols or --input"
    }
}


// ---------------------------------------------------------------------------
//  Step 1a: Bronze — fetch raw Ensembl tree, protein MSA, and CDS
// ---------------------------------------------------------------------------

process FETCH_ALIGNMENT {
    tag "${gene_symbol}"
    publishDir { "${params.outdir}/${gene_symbol.replaceAll(/[^A-Za-z0-9_]/, '_')}_${params.taxon}" },
        mode: 'copy', overwrite: true

    input:
    tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset)

    output:
    tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset), path("1_bronze"),
     emit: bronze

    script:
    def uniprot_arg = uniprot ? "--uniprot ${uniprot}" : ''
    """
    mkdir -p 1_bronze
    uv run python ${params.script_dir}/ensembl_msa_tree_extract.py \\
        --gene-symbol '${gene_symbol}' \\
        --outdir 1_bronze \\
        --taxon ${params.taxon} \\
        ${uniprot_arg}
    """
}

// ---------------------------------------------------------------------------
//  Step 1b: Bronze — Fetch 3D structure
// ---------------------------------------------------------------------------

process FETCH_STRUCTURE {
    tag "${gene_symbol}"
    publishDir { "${params.outdir}/${gene_symbol.replaceAll(/[^A-Za-z0-9_]/, '_')}_${params.taxon}" },
        mode: 'copy', overwrite: true

    input:
    tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset)

    output:
    tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset),
        path("1_bronze_struct"), emit: bronze_struct

    script:
    def uniprot_arg = uniprot ? "--uniprot ${uniprot}" : ''
    def pdb_arg     = pdb    ? "--pdb ${pdb}" : ''
    """
    mkdir -p 1_bronze_struct
    uv run python ${params.script_dir}/fetch_structure.py \\
        --gene-symbol '${gene_symbol}' \\
        --outdir 1_bronze_struct \\
        --max-pdb ${params.max_pdb} \\
        ${uniprot_arg} \\
        ${pdb_arg}
    """
}

// ---------------------------------------------------------------------------
//  Step 2b: Silver — Renumber Bronze PDBs to UniProt numbering (SIFTS)
// ---------------------------------------------------------------------------

process RENUMBER_STRUCTURES {
    tag "${gene_symbol}"
    publishDir { "${params.outdir}/${gene_symbol.replaceAll(/[^A-Za-z0-9_]/, '_')}_${params.taxon}" },
        mode: 'copy', overwrite: true

    input:
    tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset),
        path(bronze_struct)

    output:
    tuple val(gene_symbol), path("2_silver_struct"), emit: silver_struct

    script:
    def uniprot_arg = uniprot ? "--uniprot ${uniprot}" : ''
    """
    mkdir -p 2_silver_struct
    uv run python ${params.script_dir}/renumber_structures.py \\
        --gene-symbol '${gene_symbol}' \\
        --indir ${bronze_struct} \\
        --outdir 2_silver_struct \\
        ${uniprot_arg}
    """
}

// ---------------------------------------------------------------------------
//  Step 2a: Silver — transform Bronze into a CodeML-ready alignment + tree
// ---------------------------------------------------------------------------

process TRANSFORM_ALIGNMENT_SILVER {
    tag "${gene_symbol}"
    publishDir { "${params.outdir}/${gene_symbol.replaceAll(/[^A-Za-z0-9_]/, '_')}_${params.taxon}" },
    mode: 'copy'

    input:
    tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset), path(bronze)

    output:
    tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset), path("2_silver"),
    emit: silver

    script:
    def uniprot_arg = uniprot ? "--uniprot ${uniprot}" : ''
    """
    mkdir -p 2_silver
    uv run python ${params.script_dir}/ensembl_msa_tree_transform.py \\
        --gene-symbol '${gene_symbol}' \\
        --indir ${bronze} \\
        --outdir 2_silver \\
        ${uniprot_arg}
    """
}


// ---------------------------------------------------------------------------
//  Merge Step 1 outputs
// ---------------------------------------------------------------------------

// process MERGE_PREPARE {
//     tag "${gene_symbol}"
//
//     input:
//     tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset), path("silver"), path("struct_workdir")
//
//     output:
//     tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset), path("silver"), emit: workdir
//
//     script:
//     """
//     mkdir -p workdir
//     cp -a silver/* workdir/ 2>/dev/null || true
//     cp -a struct_workdir/* workdir/ 2>/dev/null || true
//     """
//}


// ---------------------------------------------------------------------------
//  Step 2b1: Silver — Run CodeML M8 (alternative model)
// ---------------------------------------------------------------------------

process RUN_CODEML_M8 {
    tag "${gene_symbol} M8"
    cpus 1
    publishDir { "${params.outdir}/${gene_symbol.replaceAll(/[^A-Za-z0-9_]/, '_')}_${params.taxon}" },
    mode: 'copy'

    input:
    tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset), path(silver)

    output:
    tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset),
    path("m8_workdir"), emit: m8_workdir

    script:
    def uniprot_arg = uniprot ? "--uniprot ${uniprot}" : ''
    """
    mkdir -p m8_workdir
    cp -a ${silver}/. m8_workdir/
    uv run python ${params.script_dir}/run_codeml_site_models.py \\
        --gene-symbol '${gene_symbol}' \\
        --outdir m8_workdir \\
        --model M8 \\
        ${uniprot_arg}
    """
}


// ---------------------------------------------------------------------------
//  Step 2b2: Silver — Run CodeML M8a (null model)
// ---------------------------------------------------------------------------

process RUN_CODEML_M8A {
    tag "${gene_symbol} M8a"
    cpus 1
    publishDir { "${params.outdir}/${gene_symbol.replaceAll(/[^A-Za-z0-9_]/, '_')}_${params.taxon}" },
    mode: 'copy'

    input:
    tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset), path(silver)

    output:
    tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset),
    path("m8a_workdir"), emit: m8a_workdir

    script:
    def uniprot_arg = uniprot ? "--uniprot ${uniprot}" : ''
    """
    mkdir -p m8a_workdir
    cp -a ${silver}/. m8a_workdir/
    uv run python ${params.script_dir}/run_codeml_site_models.py \\
        --gene-symbol '${gene_symbol}' \\
        --outdir m8a_workdir \\
        --model M8a \\
        ${uniprot_arg}
    """
}


// ---------------------------------------------------------------------------
//  Step 2b3: Extract BEB results (after both models finish)
// ---------------------------------------------------------------------------

process EXTRACT_RESULTS {
    tag "${gene_symbol}"
    publishDir { "${params.outdir}/${gene_symbol.replaceAll(/[^A-Za-z0-9_]/, '_')}_${params.taxon}" },
    mode: 'copy'

    input:
    tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset),
    path("m8_workdir"), path("m8a_workdir")

    output:
    tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset), path("2_silver"),
    emit: silver

    script:
    def uniprot_arg = uniprot ? "--uniprot ${uniprot}" : ''
    """
    mkdir -p 2_silver
    cp -a m8_workdir/* 2_silver/ 2>/dev/null || true
    cp -a m8a_workdir/* 2_silver/ 2>/dev/null || true
    uv run python ${params.script_dir}/run_codeml_site_models.py \\
        --gene-symbol '${gene_symbol}' \\
        --outdir 2_silver \\
        --extract-only \\
        ${uniprot_arg}
    """
}


// ---------------------------------------------------------------------------
//  Step 2c: Silver - Run HyPhy site models (parallel cross-check of CodeML M8)
// ---------------------------------------------------------------------------

process RUN_HYPHY {
    tag "${gene_symbol} ${params.hyphy_method}"
    publishDir { "${params.outdir}/${gene_symbol.replaceAll(/[^A-Za-z0-9_]/, '_')}_${params.taxon}" },
    mode: 'copy'
    cpus 1

    input:
    tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset), path(silver)

    output:
    tuple val(gene_symbol), path(silver), emit: silver

    script:
    def uniprot_arg = uniprot ? "--uniprot ${uniprot}" : ''
    """
    uv run python ${params.script_dir}/run_hyphy_site_models.py \\
        --gene-symbol '${gene_symbol}' \\
        --outdir 2_silver \\
        --method ${params.hyphy_method} \\
        ${uniprot_arg}
    """
}


// ---------------------------------------------------------------------------
//  Step 2d: Silver — CodeML branch-site model A along the target lineage
// ---------------------------------------------------------------------------

process RUN_CODEML_BRANCH_SITE {
    tag "${gene_symbol} branch-site"
    publishDir { "${params.outdir}/${gene_symbol.replaceAll(/[^A-Za-z0-9_]/, '_')}_${params.taxon}" },
    mode: 'copy'
    cpus 1

    input:
    tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset), path(silver)

    output:
    tuple val(gene_symbol), path(silver), emit: silver

    script:
    def uniprot_arg = uniprot ? "--uniprot ${uniprot}" : ''
    """
    uv run python ${params.script_dir}/run_codeml_branch_site_models.py \\
        --gene-symbol '${gene_symbol}' \\
        --outdir 2_silver \\
        ${uniprot_arg}
    """
}


// ---------------------------------------------------------------------------
//  Step 3: Analyse results (Jalview, PyMOL, dN/dS plot)
// ---------------------------------------------------------------------------

process ANALYZE_RESULTS {
    tag "${gene_symbol}"
    publishDir { "${params.outdir}/${gene_symbol.replaceAll(/[^A-Za-z0-9_]/, '_')}_${params.taxon}" },
    mode: 'copy'

    input:
    tuple val(gene_symbol), val(uniprot), val(pdb), val(chain), val(resi_offset),
        path("2_silver"), path("2_silver_struct")

    output:
    tuple val(gene_symbol), path("3_gold"), emit: gold

    script:
    def uniprot_arg     = uniprot ? "--uniprot ${uniprot}" : ''
    def pdb_arg         = pdb    ? "--pdb ${pdb}" : ''
    def offset_arg      = resi_offset != 0 ? "--resi-offset ${resi_offset}" : ''
    """
    mkdir -p 3_gold
    cp -a 2_silver/. 3_gold/
    cp -a 2_silver_struct/*.pdb 3_gold/ 2>/dev/null || true
    cp -a 2_silver_struct/*_structures.json 3_gold/ 2>/dev/null || true
    uv run python ${params.script_dir}/analyse_site_models.py \\
        --gene-symbol '${gene_symbol}' \\
        --outdir 3_gold \\
        --chain ${chain} \\
        --taxon ${params.taxon} \\
        --plot-threshold ${params.plot_threshold} \\
        ${uniprot_arg} \\
        ${pdb_arg} \\
        ${offset_arg}
    """
}


// ---------------------------------------------------------------------------
//  Workflow
// ---------------------------------------------------------------------------

workflow {
    ch_input = build_input_channel()

    // Step 1a: Bronze — fetch raw Ensembl tree, MSA, and CDS
    FETCH_ALIGNMENT(ch_input)

    // Step 1b: Bronze — fetch 3D structures (AlphaFold + SIFTS-ranked RCSB PDBs)
    FETCH_STRUCTURE(ch_input)

    // Step 2b (struct): Silver — renumber Bronze PDBs to UniProt numbering via SIFTS
    RENUMBER_STRUCTURES(FETCH_STRUCTURE.out.bronze_struct)

    // Step 2a: Silver — transform Bronze into a CodeML-ready alignment + tree
    TRANSFORM_ALIGNMENT_SILVER(FETCH_ALIGNMENT.out.bronze)

    ch_silver = TRANSFORM_ALIGNMENT_SILVER.out.silver

    // Step 2b: run CodeML M8 + M8a in parallel, then extract BEB
    if (params.codeml_method) {
        if (params.skip_codeml_run) {
            // Reuse silver as both m8 and m8a workdirs (expects pre-existing .mlc)
            ch_codeml_joined = ch_silver.map {
                gene_symbol, uniprot, pdb, chain, resi_offset, silver_dir ->
                tuple(gene_symbol, uniprot, pdb, chain, resi_offset, silver_dir, silver_dir)
            }
        } else {
            RUN_CODEML_M8(ch_silver)
            RUN_CODEML_M8A(ch_silver)

            // Join M8 + M8a outputs by gene_symbol, then extract
            ch_m8  = RUN_CODEML_M8.out.m8_workdir
            ch_m8a = RUN_CODEML_M8A.out.m8a_workdir.map { it -> tuple(it[0], it[5]) }

            ch_codeml_joined = ch_m8.join(ch_m8a, by: 0)
                .map { gene_symbol, uniprot, pdb, chain, resi_offset, m8_dir, m8a_dir ->
                    tuple(gene_symbol, uniprot, pdb, chain, resi_offset, m8_dir, m8a_dir)
                }
        }

        EXTRACT_RESULTS(ch_codeml_joined)

        // Step 3: analyse — join the per-gene silver_struct directory (always
        // present; may be empty) by gene_symbol.
        ch_silver_struct = RENUMBER_STRUCTURES.out.silver_struct
        ch_analyze_input = EXTRACT_RESULTS.out.silver
            .join(ch_silver_struct, by: 0)
            .map { gene_symbol, uniprot, pdb, chain, resi_offset, silver_dir, struct_dir ->
                tuple(gene_symbol, uniprot, pdb, chain, resi_offset, silver_dir, struct_dir)
            }

        ANALYZE_RESULTS(ch_analyze_input)
    }

    // Step 2c: HyPhy site models (parallel cross-check of CodeML M8)
    if (params.hyphy_method) {
        RUN_HYPHY(ch_silver)
    }

    // Step 2d: CodeML branch-site model A along the target lineage
    if (params.branch_site) {
        RUN_CODEML_BRANCH_SITE(ch_silver)
    }
}
