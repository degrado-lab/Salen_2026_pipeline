Sequence design and structure prediction (SDSP) pipeline
This pipeline implements LigandMPNN (LMPNN), Rosetta FastRelax (RSFR), RaptorX-Single (RXS), Chai1

Third party code: 
    -sequence design: LigandMPNN 
    -structure prediction: RaptorX-Single, Chai1

SDSP pipeline file system structure:
root_run_directory
    ├── input
    │   ├── 1_lmpnn_in
    │   │   └── input_pdb.pdb [STARTING STRUCTURE]
    │   ├── constraints
    │   │   └── constraints.cst [DO NOT CHANGE THIS FILENAME, CONSTRAINTS WRITTEN IN ROSETTA FORMAT]
    │   ├── omit_aa_per_residue
    │   │   └── omit_aa_per_residue.json [DO NOT CHANGE THIS FILENAME, JSON FILE FOR LMPNN]
    │   └── params
    │       └── ligand.params [THIS IS FOR ROSETTA]
    ├── output
    │   ├── 1_rs_in [LMPNN ROUND 1 OUTPUT]
    │   ├── 2_lmpnn_in [RSFR ROUND 1 OUTPUT]
    │   ├── 2_rs_in [LMPNN ROUND 2 OUTPUT]
    │   ├── 3_lmpnn_in [RSFR ROUND 2 OUTPUT]
    │   ├── 3_rsfr_in [LMPNN ROUND 3 OUTPUT]
    │   ├── analysis [ANALYSIS CSV FILES WILL BE HERE]
    │   ├── chai1_predictions
    │   │   ├── input_fastas [CHAI1 PREDICTION INPUT FASTAS]
    │   │   └── output [CHAI1 PREDICTION OUTPUT]
    │   ├── final_out [RSFR ROUND 3 OUTPUT]
    │   ├── ipynb_job_scripts [SGE SCHEDULER JOB FILES]
    │   ├── rxs_predictions [RXS PREDICTION OUTPUT]
            └── fastas [RXS PREDICTION INPUT FASTAS]
    │   └── wynton [SGE JOB OUTPUT FILES]
    └── scripts
        ├── chai1_scripts
        │   ├── helpers
        │   │   ├── analyze_chai1_helpers.py
        │   │   ├── find_fasta_file.py
        │   │   ├── prep_and_submit_chai1_helpers.py
        │   │   └── run_chai1.py
        │   ├── analyze_chai1.py
        │   └── prep_submit_analyze_chai1.py 
        └── lmpnn_rsfr_rxs_scripts
            ├── create_pdb_json.py
            ├── fast_relax_packstat_csts.xml
            ├── filter_pdbs.py
            ├── find_pdb_filename.py
            ├── lmpnn_rsfr_rxs_analysis.py
            ├── mpnn_out_to_raptor_x_input.py
            ├── mpnn_out_to_single_file_fasta.py
            ├── pdb_tag.py
            ├── remove_packed_tag.py
            └── rename_raptor_x_single.py