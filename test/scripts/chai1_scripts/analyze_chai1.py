import sys
import os
from datetime import datetime
import helpers.analyze_chai1_helpers as hp
def main(chai1_input_fasta_dir, chai1_root_output_dir, rosetta_model_dir, run_analysis_dir):
    
    """
    This chunk goes through an input directory of chai1-compatible single sequences fasta files
    and does the following for each sequence:
      - Copies corresponding rosetta model to chai1 output directory [copy_file]
      - Appends sequence information to output .cif filenames [rename_cif_files]
    """
    
    for filename in os.listdir(chai1_input_fasta_dir):
        if filename.endswith(".fasta"):
          base_filename = os.path.splitext(os.path.basename(filename))[0]
          output_dir = os.path.join(chai1_root_output_dir, base_filename)
          hp.rename_cif_and_npz_files(output_dir)
          hp.copy_contains_filename(rosetta_model_dir, output_dir, base_filename)

    directory_name = os.path.basename(os.path.normpath(chai1_root_output_dir))
    output_csv = os.path.join(run_analysis_dir, f"chai1_analysis_{directory_name}.csv")
    
    # Collect results from the analysis
    results = hp.analyze_subdirectory(chai1_root_output_dir)
    # Write results to CSV
    hp.write_results_to_csv(results, output_csv)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python main_script.py <target_designs_csv> <source_fasta_file_dir> <chai1_input_fasta_dir> <chai1_ligand_name> <chai1_ligand_smiles> <chai1_root_output_dir> <rosetta_model_dir>")
        sys.exit(1)

    chai1_input_fasta_dir = sys.argv[1] #you make this directory under input folder
    chai1_root_output_dir = sys.argv[2]
    rosetta_model_dir = sys.argv[3]
    run_analysis_dir = sys.argv[4]

    main(chai1_input_fasta_dir,chai1_root_output_dir, rosetta_model_dir, run_analysis_dir)
