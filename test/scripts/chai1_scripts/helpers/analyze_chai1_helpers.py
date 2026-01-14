import os
import shutil
import numpy as np
import subprocess
import pandas as pd
import csv
from Bio.PDB import PDBParser, MMCIFParser, Superimposer
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

def copy_contains_filename(source_directory, target_directory, filename):
    # Remove the last underscore and the following number
    base_name = filename.rsplit('_', 1)[0]
    
    # Iterate through each file in the source directory
    for file_name in os.listdir(source_directory):
        # Check if the base name is contained within the file name
        if base_name in file_name:
            source_file_path = os.path.join(source_directory, file_name)
            
            # Construct the destination path in the target directory
            destination_path = os.path.join(target_directory, file_name)
            
            # Copy the file to the target directory
            shutil.copy(source_file_path, destination_path)
            print(f"Copied {source_file_path} to {destination_path}")
            return
    
    print(f"copy_contains_filename function to get rosetta model: no matching file found for {filename} in source directory")

def rename_cif_and_npz_files(directory):
    # Get the base name of the directory
    dir_name = os.path.basename(directory)
    
    # Iterate through each file in the given directory
    for file_name in os.listdir(directory):
        # Check if the file is a .cif file
        if file_name.endswith('.cif') or file_name.endswith('.npz'):
            # Construct the new file name by appending the directory name plus an underscore
            new_file_name = f"{dir_name}_{file_name}"
            
            # Construct the full paths to the old and new file names
            old_file_path = os.path.join(directory, file_name)
            new_file_path = os.path.join(directory, new_file_name)
            
            # Rename the file
            os.rename(old_file_path, new_file_path)
            print(f"Renamed {old_file_path} to {new_file_path}")

def calculate_rmsd(pdb_file, cif_file):
    # Parse the PDB file
    pdb_parser = PDBParser(QUIET=True)
    pdb_structure = pdb_parser.get_structure('PDB_structure', pdb_file)
    
    # Parse the CIF file
    cif_parser = MMCIFParser(QUIET=True)
    cif_structure = cif_parser.get_structure('CIF_structure', cif_file)
    
    # Extract the first model and chain from each structure
    pdb_model = pdb_structure[0]
    cif_model = cif_structure[0]
    pdb_chain = pdb_model.child_list[0]
    cif_chain = cif_model.child_list[0]
    
    # Extract the atomic coordinates
    pdb_atoms = [atom for atom in pdb_chain.get_atoms() if atom.get_id() == 'CA']
    cif_atoms = [atom for atom in cif_chain.get_atoms() if atom.get_id() == 'CA']
    
    # Ensure both structures have the same number of atoms
    if len(pdb_atoms) != len(cif_atoms):
        raise ValueError("Structures have different number of atoms")
    
    # Superimpose the structures
    super_imposer = Superimposer()
    super_imposer.set_atoms(pdb_atoms, cif_atoms)
    super_imposer.apply(cif_model.get_atoms())
    
    # Calculate RMSD
    rmsd = super_imposer.rms
    
    return rmsd

def calculate_avg_plddt(cif_file):
    mmcif_dict = MMCIF2Dict(cif_file)
    
    # Extract pLDDT scores
    plddt_scores = []
    if '_ma_qa_metric_local.metric_id' in mmcif_dict and '_ma_qa_metric_local.metric_value' in mmcif_dict:
        metric_ids = mmcif_dict['_ma_qa_metric_local.metric_id']
        metric_values = mmcif_dict['_ma_qa_metric_local.metric_value']
        
        for metric_id, metric_value in zip(metric_ids, metric_values):
            if metric_id == '1':  # Assuming '1' corresponds to pLDDT
                plddt_scores.append(float(metric_value))
    avg_plddt = np.mean(plddt_scores)
    
    return avg_plddt

def get_ptm_iptm(npz_file):
    npz_scores = np.load(npz_file)
    ptm = npz_scores['ptm']
    ptm_float_value = float(ptm[0])
    iptm = npz_scores['iptm']
    iptm_float_value = float(iptm[0])
    return ptm_float_value, iptm_float_value

def get_per_chain_ptm_iptm(npz_file):
    npz_scores = np.load(npz_file)
    per_chain_ptm = npz_scores['ptm']
    per_chain_iptm = npz_scores['iptm']
    return per_chain_ptm, per_chain_iptm

def safe_mean(values):
    print(f"Calculating mean for values: {values}")
    mean_value = np.mean(values) if len(values) > 0 else float('nan')
    print(f"Mean value: {mean_value}")
    return mean_value

def safe_std(values):
    print(f"Calculating std for values: {values}")
    std_value = np.std(values) if len(values) > 0 else float('nan')
    print(f"Standard deviation: {std_value}")
    return std_value

def analyze_subdirectory(directory):
    results = []

    for root, dirs, files in os.walk(directory):
        if root == directory:
            continue
        
        if len(files) <= 1:
            print(f"Skipping directory with 1 file or less: {root}")
            continue

        pdb_file = None
        cif_files = []
        npz_files = []

        for file in files:
            if file.endswith('.pdb'):
                pdb_file = os.path.join(root, file)
            elif file.endswith('.cif'):
                cif_files.append(os.path.join(root, file))
            elif file.endswith('.npz'):
                npz_files.append(os.path.join(root, file))
        
        if pdb_file and cif_files and npz_files:
            rmsd_values = []
            plddt_values = []
            ptm_values = []
            iptm_values = []

            for cif_file in cif_files:
                try:
                    rmsd_value = calculate_rmsd(pdb_file, cif_file)
                    print(f"RMSD value for {cif_file}: {rmsd_value}")
                    rmsd_values.append(rmsd_value)
                    plddt_value = calculate_avg_plddt(cif_file)
                    print(f"pLDDT value for {cif_file}: {plddt_value}")
                    plddt_values.append(plddt_value)
                except Exception as e:
                    print(f"Error calculating RMSD or pLDDT for {cif_file}: {e}")
                    rmsd_values.append(float('nan'))
                    plddt_values.append(float('nan'))

            for npz_file in npz_files:
                try:
                    print(f"Loading .npz file: {npz_file}")
                    ptm_value, iptm_value = get_ptm_iptm(npz_file)
                    print(f"PTM value for {npz_file}: {ptm_value}, iPTM value: {iptm_value}")
                    ptm_values.append(ptm_value)
                    iptm_values.append(iptm_value)
                except Exception as e:
                    print(f"Error calculating PTM or iPTM for {npz_file}: {e}")
                    ptm_values.append(float('nan'))
                    iptm_values.append(float('nan'))

            if rmsd_values and plddt_values and ptm_values and iptm_values:

                # Ensure values are passed correctly
                print(f"Calculating mean and std for RMSD values: {rmsd_values}")
                avg_rmsd = safe_mean([value for value in rmsd_values if isinstance(value, (int, float))])
                std_rmsd = safe_std([value for value in rmsd_values if isinstance(value, (int, float))])
                
                print(f"Calculating mean and std for pLDDT values: {plddt_values}")
                avg_plddt = safe_mean([value for value in plddt_values if isinstance(value, (int, float))])
                std_plddt = safe_std([value for value in plddt_values if isinstance(value, (int, float))])
                
                print(f"Calculating mean and std for PTM values: {ptm_values}")
                avg_ptm = safe_mean([value for value in ptm_values if isinstance(value, (int, float))])
                std_ptm = safe_std([value for value in ptm_values if isinstance(value, (int, float))])
                
                print(f"Calculating mean and std for iPTM values: {iptm_values}")
                avg_iptm = safe_mean([value for value in iptm_values if isinstance(value, (int, float))])
                std_iptm = safe_std([value for value in iptm_values if isinstance(value, (int, float))])
                
                results.append({
                    'Subdirectory': os.path.basename(root),
                    'Average RMSD': avg_rmsd,
                    'RMSD Stdev': std_rmsd,
                    'Average pLDDT': avg_plddt,
                    'pLDDT Stdev': std_plddt,
                    'Average PTM': avg_ptm,
                    'PTM Stdev': std_ptm, 
                    'Average iPTM': avg_iptm,
                    'iPTM Stdev': std_iptm, 
                })
            else:
                print("One or more arrays are empty, skipping calculations.")
    return results

def write_results_to_csv(results, output_csv):
    # Sort results by 'Average RMSD' from lowest to highest
    sorted_results = sorted(results, key=lambda x: x['Average RMSD'])
    
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['Subdirectory', 'Average RMSD', 'RMSD Stdev', 'Average pLDDT', 'pLDDT Stdev',
                      'Average PTM', 'PTM Stdev', 'Average iPTM', 'iPTM Stdev']
        
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        
        writer.writeheader()
        for result in sorted_results:
            # Round the values to 2 decimal places
            result['Average RMSD'] = round(result['Average RMSD'], 2)
            result['RMSD Stdev'] = round(result['RMSD Stdev'], 2)
            result['Average pLDDT'] = round(result['Average pLDDT'], 2)
            result['pLDDT Stdev'] = round(result['pLDDT Stdev'], 2)
            result['Average PTM'] = round(result['Average PTM'], 2)
            result['PTM Stdev'] = round(result['PTM Stdev'], 2)
            result['Average iPTM'] = round(result['Average iPTM'], 2)
            result['iPTM Stdev'] = round(result['iPTM Stdev'], 2)
            writer.writerow(result)