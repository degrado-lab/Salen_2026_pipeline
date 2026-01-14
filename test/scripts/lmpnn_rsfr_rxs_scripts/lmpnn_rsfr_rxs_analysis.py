'''
rmsd (x)
side chain rmsd
avg plddt (x)
netcharge (x)
ligand sasa (x)
rosetta score/packstat (ip)
seq alignment/logo 
'''
import os
import sys
import glob
import re
from Bio import PDB
import pandas as pd
import tempfile
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

#this is a helper function that extracts the fasta sequence of a given design
def parse_fasta(fasta_file, string_variable):
  """
  Parses a FASTA file and returns a dictionary of sequences where the headers
  contain the specified string variable.

  Args:
    fasta_file: The path to the FASTA file.
    string_variable: The string variable to search for in the FASTA headers.

  Returns:
    A dictionary of sequences where the keys are the headers containing the string
    variable and the values are the corresponding sequences.
  """
  sequences = {}
  current_header = None
  with open(fasta_file) as f:
    for line in f:
      line = line.rstrip()
      if line.startswith('>'):
        current_header = line.strip()[1:]
        sequences[current_header] = ""
      elif current_header:
        sequences[current_header] += line.strip()
  return {k: v for k, v in sequences.items() if string_variable in k}

#calculates RMSD between rosetta model and omegafold prediction
def calculate_rmsd(pdb_file1, pdb_file2):
    try:
        # Parse PDB files
        parser = PDB.PDBParser(QUIET=True)
        structure1 = parser.get_structure("structure1", pdb_file1)
        structure2 = parser.get_structure("structure2", pdb_file2)

        # Extract Cα atoms from the structures
        atoms1 = [atom for atom in structure1.get_atoms() if atom.get_name() == "CA"]
        atoms2 = [atom for atom in structure2.get_atoms() if atom.get_name() == "CA"]

        # Check if the number of Cα atoms in both structures is the same
        if len(atoms1) != len(atoms2):
            raise ValueError("The number of Cα atoms in the two structures is different.")

        # Superimpose structures based on Cα atoms
        superimposer = PDB.Superimposer()
        superimposer.set_atoms(atoms1, atoms2)

        # Calculate the RMSD value
        rmsd = superimposer.rms

        return rmsd
    except Exception as e:
        # Log the mismatched comparison to a separate file
        with open("mismatched_comparisons.txt", "a") as log_file:
            log_file.write(f"Error comparing {pdb_file1} and {pdb_file2}: {str(e)}\n")
        return None

#calculates average plddt of each omegafold model 
def calculate_average_plddt(pdb_file):
    # Extract pLDDT values from the B-factor column
    plddt_values = []
    with open(pdb_file, "r") as pdb_handle:
        for line in pdb_handle:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                b_factor = float(line[60:66])
                plddt_values.append(b_factor)

    # Calculate the average pLDDT
    if plddt_values:
        average_plddt = sum(plddt_values) / len(plddt_values)
    else:
        average_plddt = 0

    return average_plddt

def calculate_net_charge(protein_seq, seq_fasta_file):
   sequence_and_id = parse_fasta(seq_fasta_file, protein_seq)
   # print(sequence_and_id) troubleshooting print statement
   sequence = sequence_and_id[protein_seq][0:]
   #print(sequence) troubleshooting print statement

    # count number of positively charged residues (K and R)
   pos_count = len(re.findall('[KR]', sequence))
    # count number of negatively charged residues (D and E) 
   neg_count = len(re.findall('[DE]', sequence))
    # calculate net charge
   net_charge = pos_count - neg_count
   return(net_charge)

def calculate_delta_sasa_binding_affinity(pdb_file_path):

#this chunk makes the temporary pdb file of just the free ligand
    hetatm_lines = []
    conect_lines = []

    with open(pdb_file_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('HETATM') or line.startswith('CONECT'):
                if line.startswith('HETATM'):
                    hetatm_lines.append(line)
                elif line.startswith('CONECT'):
                    conect_lines.append(line)

    # Create a temporary file
    temp_file_path = tempfile.mktemp(suffix=".pdb", prefix="temp_", dir=".")

    # Write HETATM and CONECT lines to temporary file
    with open(temp_file_path, 'w') as temp_file:
        temp_file.writelines(hetatm_lines)
        temp_file.writelines(conect_lines)

#sasa calculation set-up
    probe_radius = 1.4  # Modify this value as needed

    # Create a PDB parser
    p = PDBParser(QUIET=1)

    # Initialize the SASA calculator with the defined probe radius
    sr = ShrakeRupley(probe_radius=probe_radius)

#free ligand sasa calculation
    # Load the free ligand PDB structure
    free_ligand_structure = p.get_structure(os.path.basename(temp_file_path), temp_file_path)
    sr.compute(free_ligand_structure, level="A")
    
    free_ligand_sasa = sum(atom.sasa for atom in free_ligand_structure.get_atoms())

#ligand in complex sasa calculation
    # Load the ligand in complex PDB structure
    ligand_in_complex_structure = p.get_structure(os.path.basename(pdb_file_path), pdb_file_path)
    sr.compute(ligand_in_complex_structure, level="A")

    # Identify ligand residue names (HETATOMs)
    ligand_residue_name = None

    with open(pdb_file_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('HETATM'):
                # Find the residue name in the line (columns 17-20, or any other appropriate location)
                ligand_residue_name = line[17:20].strip()  # Adjust these indices based on your actual PDB format
                break

    if ligand_residue_name is None:
        raise ValueError("No ligand residue found in the PDB file.")

    # Filter atoms belonging to the ligand residue
    ligand_atoms = [atom for atom in ligand_in_complex_structure.get_atoms() if atom.get_parent().resname == ligand_residue_name]

    # Collect SASA values for ligand atoms and calculate the total SASA of ligand
    ligand_in_complex_sasa = 0
    for atom in ligand_atoms:
        ligand_in_complex_sasa += atom.sasa
    
    #delta sasa calculations   
    delta_sasa = free_ligand_sasa - ligand_in_complex_sasa

    
    #Calculate Kd based on lig_eff = -delta_g/delta sasa
    # following ligand efficiency constants are from this paper: https://doi.org/10.1021/ci200612f
    #units are kcal/mol x Angstroms^2
    # lig_eff_hard = -0.117 # absolute limit of ligand efficiency
    #lig_eff_soft = -0.051 # 95% of paper's data set is within this limit
    lig_eff_avg = -0.025 # average ligand efficiency of dataset 
    
    # following ligand efficiency constant from this paper: https://doi.org/10.1002/anie.200200565
    #lig_eff_houk = -0.007 

    delta_g = -lig_eff_avg * delta_sasa

    ideal_gas_constant = .001987  # kcal/mol*K
    room_temperature = 298  # Kelvin
    kd = np.exp(-delta_g / (ideal_gas_constant * room_temperature)) * 1e9 #nM

    os.remove(temp_file_path)
    return free_ligand_sasa, ligand_in_complex_sasa, delta_sasa

def extract_rosetta_scores(rosetta_output_file, design_name):
     # Open the text file in read mode
  with open(rosetta_output_file, "r") as text_file:
    # Read the entire text data
    text_data = text_file.read()

  # Use the previous logic to find and return the value
  return find_string_value_internal(text_data, design_name)

# Internal function to avoid code duplication for extract_rosetta_scores
def find_string_value_internal(text_data, design_name):
  lines = text_data.splitlines()

  for line in lines:
    columns = line.split("\t")

    if design_name in columns:
      return columns[0], columns[1]

  return None

def main(models_dir, predictions_dir):
    '''
    the following chunk is setting up a pandas dataframe to store the analysis data
    and defining global variables for the script 
    ''' 
    df = pd.DataFrame()
    # Extract the prediction directory name to use in the output filenames
    model_files = [file for file in os.listdir(models_dir) if file.endswith(".pdb")]
    total_files = len(model_files)
    print(f"Total files to compare: {total_files}")
    # Find fasta file in rosetta directory 
    find_fasta = glob.glob(os.path.join(models_dir, "*.fasta"))
    chosen_fasta_file = find_fasta[0]
    # Find total_score.txt file in rosetta directory 
    find_score = glob.glob(os.path.join(models_dir, "*.txt"))
    score_file = find_score[0]

    for i, model_filename in enumerate(model_files, start=1):
        try:
            parts = model_filename.rsplit("_", 1)
            print (parts)
            predictions_filename = parts[0] + '.pdb'
            print (f'prediction filename: {predictions_filename}')
            predictions_filepath = os.path.join(predictions_dir, predictions_filename)
            model_filepath = os.path.join(models_dir, model_filename) 
            model_design_name = model_filename[:-4]
            prediction_design_name = parts[0]

            if os.path.isfile(model_filepath):
                # rmsd
                rmsd_value = calculate_rmsd(model_filepath, predictions_filepath)
                # plddt
                average_plddt = calculate_average_plddt(predictions_filepath)
                #net charge
                net_charge = calculate_net_charge(prediction_design_name, chosen_fasta_file)
                #sasa
                free_ligand_sasa, ligand_in_complex_sasa, delta_sasa = calculate_delta_sasa_binding_affinity(model_filepath)
                #rosetta scores
                rosetta_total_score, packstat = extract_rosetta_scores(score_file, model_design_name)

                row = {}
                row['Design'] = model_design_name
                row['pLDDT'] = average_plddt
                row['RMSD'] = rmsd_value
                row['Rosetta total score'] = rosetta_total_score
                row['Rosetta packstat'] = packstat
                row['Net charge'] = net_charge
                row['FL SASA'] = free_ligand_sasa
                row['LIC SASA'] = ligand_in_complex_sasa
                row['DeltaSASA'] = delta_sasa

                # Concatenate the row to the DataFrame
                df = pd.concat([df, pd.DataFrame.from_dict([row])], ignore_index=True)

                # Print progress
                print(f"Processed file {i}/{total_files}: {model_filename}")
            else:
                print(f"Model file {model_filename} not found.")
        except Exception as e:
            print(f"An error occurred while processing {model_filename}: {e}")
            # Optionally, add logging or other error handling here

    # Round all numeric columns to two decimal places
    df_rounded = df.round(2)

    # Export the rounded DataFrame to CSV
    df_rounded.to_csv('lmpnn_rsfr_rxs_analysis.csv', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py rosetta_models_directory omegafold_predictions_directory")
        sys.exit(1)

    rosetta_models_directory = sys.argv[1]
    omegafold_predictions_directory = sys.argv[2]

    if not os.path.isdir(rosetta_models_directory) or not os.path.isdir(omegafold_predictions_directory):
        print("Invalid directory paths provided.")
        sys.exit(1)

    main(rosetta_models_directory, omegafold_predictions_directory)
