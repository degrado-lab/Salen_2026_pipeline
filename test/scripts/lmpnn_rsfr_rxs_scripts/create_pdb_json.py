import os
import sys
import json

def create_pdb_json(pdb_directory):

    # Create an empty dictionary to store file paths
    pdb_paths = {}

    # Iterate through each file in the directory
    for filename in os.listdir(pdb_directory):
        if filename.endswith(".pdb"):
            filepath = os.path.join(pdb_directory, filename)
            pdb_paths[filepath] = ""  # Store the path with an empty value

    # Create the JSON file in the same directory as the PDB files
    with open(os.path.join(pdb_directory, "pdb_file_paths.json"), "w") as json_file:
        json.dump(pdb_paths, json_file, indent=4)
    print("JSON file generated successfully in the PDB directory!")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python create_pdb_json.py <pdb_directory>")
    sys.exit(1)

    input_directory = sys.argv[1]
    create_pdb_json(input_directory)
