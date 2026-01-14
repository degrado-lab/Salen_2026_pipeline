import os
import sys

def rename_pdb_files(directory):
  """Renames PDB files in a directory by removing the last underscore and everything after it.

  Args:
    directory: The directory containing the PDB files.
  """

  for file in os.listdir(directory):
    if file.endswith(".pdb"):
      new_name = file[:file.rfind("_")] + ".pdb"
      old_path = os.path.join(directory, file)
      new_path = os.path.join(directory, new_name)
      os.rename(old_path, new_path)
      print(f"Renamed {file} to {new_name}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(1)

    input_directory = sys.argv[1]
    rename_pdb_files(input_directory)