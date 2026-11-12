
# =============================================================================
# Protein-Protein Affinity Optimisation Pipeline
# =============================================================================

# Configuration variables
PDB_ID := 8GZ5
CHAIN_ID := B
OTHER_MOL := A
DDG_THRESHOLD := -1.3
PARALLEL_JOBS := 7

# Script directory
SCRIPTS_DIR := ./scripts
PYTHON := python

# Output files
WORKDIR := ${PDB_ID}_workdir
COMMAND_LIST_FILE := ${WORKDIR}/command_list.txt
COMMAND_LIST_PAIRS_FILE := ${WORKDIR}/command_list_pairs.txt
COMMAND_LIST_TRIPLETS_FILE := ${WORKDIR}/command_list_triplets.txt

# =============================================================================
# Phony targets
# =============================================================================
.PHONY: all help clean download_pdb repair_pdb identify_interface extract_interface \
        prepare_pssm_commands run_foldx_pssm analyse_pssm_output \
        prepare_multimutant_commands run_foldx_pairs run_foldx_triplets \
        analyse_multimutant_output pssm_analysis multimutant_analysis

# =============================================================================
# Default target
# =============================================================================
.DEFAULT_GOAL := help

# =============================================================================
# Help target
# =============================================================================
help: ## Show this help message
	@echo "Protein-Protein Affinity Optimization Pipeline"
	@echo "=============================================="
	@echo ""
	@echo "Configuration:"
	@echo "  PDB_ID=$(PDB_ID)"
	@echo "  CHAIN_ID=$(CHAIN_ID)"
	@echo "  OTHER_MOL=$(OTHER_MOL)"
	@echo "  DDG_THRESHOLD=$(DDG_THRESHOLD)"
	@echo "  PARALLEL_JOBS=$(PARALLEL_JOBS)"
	@echo ""
	@echo "Available targets:"
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "  %-25s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

# =============================================================================
# Main pipeline targets
# =============================================================================
all: ## Run complete pipeline
all: download_pdb repair_pdb identify_interface extract_interface pssm_analysis multimutant_analysis

pssm_analysis: ## Run PSSM analysis workflow
pssm_analysis: prepare_pssm_commands run_foldx_pssm analyse_pssm_output

multimutant_analysis: ## Run multimutant analysis workflow
multimutant_analysis: prepare_multimutant_commands run_foldx_pairs run_foldx_triplets analyse_multimutant_output

# =============================================================================
# Individual pipeline steps
# =============================================================================
download_pdb: ## Download PDB structure
	$(PYTHON) $(SCRIPTS_DIR)/download_pdb.py $(PDB_ID) $(WORKDIR)

repair_pdb:  ## Repair PDB structure using FoldX
	$(PYTHON) $(SCRIPTS_DIR)/run_foldx_repair.py $(PDB_ID).pdb $(WORKDIR)

identify_interface:  ## Identify protein-protein interface
	$(PYTHON) $(SCRIPTS_DIR)/run_foldx_analysecomplex.py $(PDB_ID)_Repair.pdb $(OTHER_MOL) $(CHAIN_ID) $(WORKDIR)

extract_interface: ## Extract interface residues
	$(PYTHON) $(SCRIPTS_DIR)/extract_interface.py $(PDB_ID) $(CHAIN_ID) $(WORKDIR)

prepare_pssm_commands: ## Generate PSSM scanning commands
	$(PYTHON) $(SCRIPTS_DIR)/generate_pssm_commands_file.py $(PDB_ID) $(OTHER_MOL) $(CHAIN_ID) $(WORKDIR)

run_foldx_pssm: ## Run FoldX PSSM scanning in parallel
	@if [ -f $(COMMAND_LIST_FILE) ]; then \
		cat $(COMMAND_LIST_FILE) | xargs -P $(PARALLEL_JOBS) -I {} sh -c '{}'; \
	else \
		echo "Error: $(COMMAND_LIST_FILE) not found. Run 'make prepare_pssm_commands' first."; \
		exit 1; \
	fi

analyse_pssm_output: ## Parse and plot PSSM results
	$(PYTHON) $(SCRIPTS_DIR)/parse_pssm_output.py $(PDB_ID) $(CHAIN_ID) $(WORKDIR)
	$(PYTHON) $(SCRIPTS_DIR)/plot_pssm_heatmap.py $(PDB_ID) $(WORKDIR)

prepare_multimutant_commands: ## Generate multimutant commands
	$(PYTHON) $(SCRIPTS_DIR)/generate_multimutant_commands_file.py $(PDB_ID) $(DDG_THRESHOLD) $(WORKDIR)

run_foldx_pairs: ## Run FoldX for pairs in parallel
	@if [ -f $(COMMAND_LIST_PAIRS_FILE) ]; then \
		cat $(COMMAND_LIST_PAIRS_FILE) | xargs -P $(PARALLEL_JOBS) -I {} sh -c '{}'; \
	else \
		echo "Error: $(COMMAND_LIST_PAIRS_FILE) not found. Run 'make prepare_multimutant_commands' first."; \
		exit 1; \
	fi

run_foldx_triplets: ## Run FoldX for triplets in parallel
	@if [ -f $(COMMAND_LIST_TRIPLETS_FILE) ]; then \
		cat $(COMMAND_LIST_TRIPLETS_FILE) | xargs -P $(PARALLEL_JOBS) -I {} sh -c '{}'; \
	else \
		echo "Error: $(COMMAND_LIST_TRIPLETS_FILE) not found. Run 'make prepare_multimutant_commands' first."; \
		exit 1; \
	fi

analyse_multimutant_output:  ## Analyse multimutant results
	$(PYTHON) $(SCRIPTS_DIR)/analyse_multivariants.py $(PDB_ID) $(WORKDIR)

# =============================================================================
# Utility targets
# =============================================================================
clean: ## Clean generated files
	@echo "Cleaning generated files..."
	@rm -f $(COMMAND_LIST_FILE) $(COMMAND_LIST_PAIRS_FILE) $(COMMAND_LIST_TRIPLETS_FILE)
	@rm -rf output_* *.log *.tmp
	@echo "Clean completed."

variables: ## Show current variable values
	@echo "Current configuration:"
	@echo "  PDB_ID=$(PDB_ID)"
	@echo "  CHAIN_ID=$(CHAIN_ID)"
	@echo "  OTHER_MOL=$(OTHER_MOL)"
	@echo "  DDG_THRESHOLD=$(DDG_THRESHOLD)"
	@echo "  PARALLEL_JOBS=$(PARALLEL_JOBS)"
	@echo "  SCRIPTS_DIR=$(SCRIPTS_DIR)"
	@echo "  PYTHON=$(PYTHON)"
