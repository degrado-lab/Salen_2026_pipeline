'''
this script preps input fasta files and output directories for chai1 run
'''
import sys
import os
import helpers.prep_and_submit_chai1_helpers as hp

def prep_chai1(csv_of_target_fastas, source_fasta_file_dir, chai1_input_fasta_dir, chai1_run_root_output_dir, 
         chai1_ligand_name, chai1_ligand_smiles):
    
    """
    This function creates chai1 compatible fasta files by:
      - Finds and copies single sequence fasta from 
      lmpnn_rs pipeline structure_predictions/fasta/ [find_file_with_string]
      - Modifies protein sequence header line and appends [modify_fasta_file]
    """
    hp.prep_holo_fasta_files(csv_of_target_fastas, source_fasta_file_dir, chai1_input_fasta_dir,
                   chai1_ligand_name, chai1_ligand_smiles)
    
    hp.prep_apo_fasta_files(chai1_input_fasta_dir)
    
    """
    This function tags fasta files with a _# barcode so it can be recognized in array format job submission
    """
    hp.tag_fastas(chai1_input_fasta_dir)

    """
    creates output directory for each fasta sequence in the chai1_run_root_output_dir
    """
    for filename in os.listdir(chai1_input_fasta_dir):
      if filename.endswith(".fasta"):
        base_filename = os.path.splitext(os.path.basename(filename))[0]
        output_dir = os.path.join(chai1_run_root_output_dir, base_filename)
        os.makedirs(output_dir, exist_ok=True)
        print(f"{filename} directory has been made")

def submit_chai1(chai1_run_input_fasta_dir, chai1_run_root_output_dir, wynton_output_dir, ipynb_job_dir, chai1_scripts_path):
    number_of_runs = hp.count_files_in_directory(chai1_run_input_fasta_dir)    
    ipynb_job_script = f"{ipynb_job_dir}/submit_chai1_run.sh"
    bash_script = f'''
#!/bin/bash
#$ -S /bin/bash
#$ -o {wynton_output_dir}
#$ -cwd
#$ -j y
#$ -l gpu_mem=32G
#$ -l h_rt=00:12:00 #change as needed
#$ -q gpu.q
#$ -t 1-{number_of_runs}

export CUDA_VISIBLE_DEVICES=$SGE_GPU
TRANSFORMERS_OFFLINE=1

conda activate Chai1
fasta_tag=$SGE_TASK_ID  # Get the PDB tag from the task ID
input_fasta=$(python3 {chai1_scripts_path}/helpers/find_fasta_file.py {chai1_run_input_fasta_dir}/ _$fasta_tag.fasta)
fasta_basename=$(basename "$input_fasta")
fasta_basename_no_ext="${{fasta_basename%.*}}"
output_dir={chai1_run_root_output_dir}/$fasta_basename_no_ext

python3 {chai1_scripts_path}/helpers/run_chai1.py $input_fasta $output_dir
    '''

    with open(ipynb_job_script, 'w') as file:
            file.write(bash_script)

    job_id_submit_chai1 = hp.queue_on_wynton(ipynb_job_script, wait_for_jobid=None)
    return job_id_submit_chai1

def analyze_chai1(chai1_run_input_fasta_dir, chai1_run_root_output_dir, wynton_output_dir, ipynb_job_dir, chai1_scripts_path, rosetta_model_dir, job_id_submit_chai1, run_analysis_dir):
    ipynb_job_script = f"{ipynb_job_dir}/analyze_chai1_run.sh"
    bash_script = f'''
#!/bin/bash
#$ -S /bin/bash
#$ -o {wynton_output_dir}
#$ -cwd
#$ -j y
#$ -l h_rt=00:20:00 #change as needed
#$ -t 1

export CUDA_VISIBLE_DEVICES=$SGE_GPU
TRANSFORMERS_OFFLINE=1

conda activate chai1_analysis
fasta_tag=$SGE_TASK_ID  # Get the PDB tag from the task ID
input_fasta=$(python3 {chai1_scripts_path}/helpers/find_fasta_file.py {chai1_run_input_fasta_dir}/ _$fasta_tag.fasta)
fasta_basename=$(basename "$input_fasta")
fasta_basename_no_ext="${{fasta_basename%.*}}"
output_dir={chai1_run_root_output_dir}/$fasta_basename_no_ext

python3 {chai1_scripts_path}/analyze_chai1.py {chai1_input_fasta_dir} {chai1_root_output_dir} {rosetta_model_dir} {run_analysis_dir}
    '''

    with open(ipynb_job_script, 'w') as file:
            file.write(bash_script)

    job_id_analyze_chai1 = hp.queue_on_wynton(ipynb_job_script, wait_for_jobid=job_id_submit_chai1)


if __name__ == "__main__":
    if len(sys.argv) != 12:
        print("Usage: python prep_chai1.py <target_designs_csv> <source_fasta_file_dir> <chai1_input_fasta_dir> <chai1_ligand_name> <chai1_ligand_smiles>")
        sys.exit(1)

    target_designs_csv = sys.argv[1] #designs you want chai1 to predict
    source_fasta_file_dir = sys.argv[2] #fastas directory under structure prediction
    chai1_input_fasta_dir = sys.argv[3] #input folder containing input fasta files
    chai1_root_output_dir = sys.argv[4] #root output directory where design specific output directories are
    chai1_ligand_name = sys.argv[5] #anything but typically use PDB 3 character ID
    chai1_ligand_smiles = sys.argv[6] #ligand smiles string, typically want to use canonical smiles string
    ipynb_job_dir = sys.argv[7]
    chai1_scripts_path = sys.argv[8]
    wynton_output_dir = sys.argv[9]
    rosetta_model_dir = sys.argv[10]
    run_analysis_dir = sys.argv[11]

    prep_chai1(target_designs_csv, source_fasta_file_dir, chai1_input_fasta_dir, chai1_root_output_dir,
         chai1_ligand_name, chai1_ligand_smiles)

    job_id_submit_chai1 = submit_chai1(chai1_input_fasta_dir, chai1_root_output_dir, wynton_output_dir, ipynb_job_dir, chai1_scripts_path)


    analyze_chai1(chai1_input_fasta_dir, chai1_root_output_dir, wynton_output_dir, ipynb_job_dir, chai1_scripts_path, rosetta_model_dir, job_id_submit_chai1, run_analysis_dir)