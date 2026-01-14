import os
import csv
import shutil
import subprocess

def find_file_with_string(directory, search_string):
  """
  Finds the first file in the given directory whose filename contains the specified string.

  Args:
      directory: The directory path to search in.
      search_string: The string to search for in the filenames.

  Returns:
      The path to the first matching file, or None if no file is found.
""" 
  for filename in os.listdir(directory):
      if search_string in filename:  # Check if the filename contains the string
          filepath = os.path.join(directory, filename)
          if os.path.isfile(filepath):  # Ensure it's a file, not a directory
              print(f"successful_file:{filepath}")
              return filepath

  return None  # No matching file found

def chai1_modify_fasta_file(file_path, ligand_name, smiles_string):
    """
    Modifies a file by adding a header and a ligand line.

    Args:
    file_path: Path to the file to be modified.
    smiles_string: The SMILES string to be added to the ligand line.

    Returns:
    None
    """
    with open(file_path, 'r+') as f:
        # Read the first two lines
        line1 = f.readline()
        line2 = f.readline()

        # Modify the first line
        first_line = line1.strip('>')
        modified_line1 = f">protein|name={first_line.strip()}\n"

        # Write the modified lines and the ligand information
        f.seek(0)  # Move back to the beginning
        f.write(modified_line1)
        f.write(line2)
        f.write(f">ligand|name={ligand_name}\n")
        f.write(f"{smiles_string}\n")
        print(f"Copied fasta file has been modified")

def prep_holo_fasta_files(csv_file, structure_pred_dir, output_dir, ligand_name, ligand_smiles):
    """
    Extract fasta files based on given input CSV to output directory, limited to 100 rows after the header row.

    Args:
        csv_file: The path to the CSV file containing design names.
        structure_pred_dir: Path to single file FASTAs in structure_predictions directory.
        output_dir: The directory to copy the extracted fasta files.
        ligand_name: Name of the ligand to be used for modification.
        ligand_smiles: SMILES string of the ligand for modification.
    """
    # Open the CSV file
    with open(csv_file, 'r') as csvfile:
        # Create a CSV reader object
        reader = csv.reader(csvfile)
        # Skip the header row
        next(reader, None)  # Uncomment this line if your CSV has a header row

        # Initialize a counter for the rows processed
        row_count = 0

        # Iterate through each row, stopping after 100 rows
        for row in reader:
            if row_count >= 100:
                break

            # Access the design name
            first_item = row[0]  # Assuming the design name is in the first column
            value = first_item.rsplit('_', 1)[0] + '.fasta'
            prediction_path = find_file_with_string(structure_pred_dir, value)

            # Copy the files with preserved permissions
            try:
                shutil.copy2(prediction_path, output_dir)
                print(f"Copied fasta file: {prediction_path}")
            except Exception as e:
                print(f"Error copying files for design {value}: {e}")

            copied_file_path = os.path.join(output_dir, value)
            chai1_modify_fasta_file(copied_file_path, ligand_name, ligand_smiles)

            # Increment the row counter
            row_count += 1

def prep_apo_fasta_files(output_dir):
    """
    Copies FASTA files from an input directory, renames them, removes the last two lines,
    and inserts "_apo" at the end of ">protein|name=".

    Args:
        input_dir (str): The path to the directory containing the FASTA files.
    """
    for filename in os.listdir(output_dir):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            input_filepath = os.path.join(output_dir, filename)
            base_filename, file_extension = os.path.splitext(filename)
            new_filename = base_filename + "_apo" + file_extension
            output_filepath = os.path.join(output_dir, new_filename)

            try:
                # Copy the file
                shutil.copy2(input_filepath, output_filepath)

                # Remove the last two lines from the *copied* file and modify the protein name
                with open(output_filepath, "r+") as f:
                    lines = f.readlines()
                    if len(lines) > 2:
                        # Modify the first line if it starts with ">protein|name="
                        if lines[0].startswith(">protein|name="):
                            parts = lines[0].split("=")
                            if len(parts) > 1:
                                parts[1] = parts[1].strip() + "_apo\n"
                                lines[0] = "=".join(parts)

                        f.seek(0)
                        f.writelines(lines[:-2])
                        f.truncate()
                        print(f"Processed: {filename} -> {new_filename}")
                    else:
                        print(f"Skipped: {filename} (less than 2 lines)")

            except Exception as e:
                print(f"Error processing {filename}: {e}")

def tag_fastas(fasta_directory):
    # Initialize a counter for the number tag
    number_tag = 1

    # Iterate through each file in the directory
    for filename in os.listdir(fasta_directory):
        if filename.endswith(".fasta"):
            filepath = os.path.join(fasta_directory, filename)

            # Split filename into base and extension
            base, extension = os.path.splitext(filename)

            # Construct new filename with base, new tag, and extension
            new_filename = f"{base}_{number_tag}{extension}"
            new_filepath = os.path.join(fasta_directory, new_filename)

            # Rename the file
            os.rename(filepath, new_filepath)
            print(new_filepath)

            # Increment the tag for the next file
            number_tag += 1

    print("Fasta files successfully tagged!")

def count_files_in_directory(directory):
    """
    Counts the number of files in the specified directory.

    Parameters:
    directory (str): The path to the directory.

    Returns:
    int: The number of files in the directory.
    """
    try:
        # List all items in the directory and filter out the files
        file_count = len([name for name in os.listdir(directory) if os.path.isfile(os.path.join(directory, name))])
        return file_count
    except Exception as e:
        print(f"An error occurred: {e}")
        return 0

def queue_on_wynton(run_script, wait_for_jobid=None):
    '''
    Arguments:
        run_script: (string) 
            The path of the script to run.
        wait_for_jobid: (None or String, default 'None')
            Job ID of a previous job.
            If specified, the current job will be queued, but wait for the specified job to be completed before running.
    Returns:
        jobid: (string) 
            The job id of the queued job.
            Use this to check the status of the job using 'qstat'.
    '''
    
    #Finally, let's queue the job:
    if wait_for_jobid is None:
        command = ['qsub', '-cwd', str(run_script)] #subprocess takes commands as a list of strings
    else:
        command = ['qsub', '-cwd', '-hold_jid', str(wait_for_jobid), str(run_script)] 
    
    #execute:
    result = subprocess.run(command, capture_output=True)

    if result.returncode != 0:
        print("Error: Job not submitted. Please check the error message below.")
        error = result.stderr.strip()
        print(error.decode())
        return None
    
    output = result.stdout.strip()
    print(output.decode()) #turn from bytes to string and print
    print("Use the command 'qstat' in the console to check the status of your job.")
    jobid = output.decode().split()[2] #get the job id
    jobid = jobid.split('.')[0] #remove the suffix denoting an array job.
    return jobid