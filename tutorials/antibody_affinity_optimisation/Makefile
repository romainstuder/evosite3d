
PDB_ID = 9FWW
CHAIN_ID = B
OTHER_MOL = A
DDG_THRESHOLD = -1.5


default: all

all: download_pdb \
     repair_pdb \
	 identify_interface \
	 extract_interface \
     prepare_pssm_commands \
     run_foldx_pssm \
     analyse_pssm_output \
     prepare_multimutant_commands \
     run_foldx_pairs \
     run_foldx_triplets \
     analyse_multimutant_output \
     get_final_score

download_pdb:
	python ./scripts/download_pdb.py $(PDB_ID)

repair_pdb:
	python ./scripts/run_foldx_repair.py $(PDB_ID)

identify_interface:
	python ./scripts/run_foldx_analysecomplex.py $(PDB_ID) $(OTHER_MOL) $(CHAIN_ID)

extract_interface:
	python ./scripts/extract_interface.py $(PDB_ID) $(CHAIN_ID)

prepare_pssm_commands:
	python ./scripts/generate_pssm_commands_file.py $(PDB_ID) $(OTHER_MOL) $(CHAIN_ID)

run_foldx_pssm:
	cat command_list.txt | xargs -P 7 -I {} sh -c '{}'

analyse_pssm_output:
	python ./scripts/parse_pssm_output.py $(PDB_ID) $(CHAIN_ID)
	python ./scripts/plot_pssm_heatmap.py $(PDB_ID)

prepare_multimutant_commands:
	python ./scripts/generate_multimutant_commands_file.py $(PDB_ID) $(DDG_THRESHOLD)

run_foldx_pairs:
	cat command_list_pairs.txt | xargs -P 7 -I {} sh -c '{}'

run_foldx_triplets:
	cat command_list_triplets.txt | xargs -P 7 -I {} sh -c '{}'

analyse_multimutant_output:
	python ./scripts/analyse_multivariants.py $(PDB_ID)

get_final_score:
	python ./scripts/compute_scoring.py $(PDB_ID) $(CHAIN_ID)
