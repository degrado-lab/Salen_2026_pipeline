import os
import sys

def pdb_tag(pdb_directory):
# Initialize a counter for the number tag
    number_tag = 0

    # Iterate through each file in the directory
    for filename in os.listdir(pdb_directory):
        if filename.endswith(".pdb"):
            filepath = os.path.join(pdb_directory, filename)

            # Split filename into base and existing tag (if any)
            base, existing_tag = os.path.splitext(filename)[0].rsplit("_", 1)

            # Increment tag only if a numerical tag already exists
            if existing_tag.isdigit():
                number_tag += 1
            else:
                existing_tag = ""

            # Construct new filename with base, existing tag, new tag, and extension
            new_filename = base + "_" + existing_tag + "_" + str(number_tag) + ".pdb"
            new_filepath = os.path.join(pdb_directory, new_filename)

            # Rename the file
            os.rename(filepath, new_filepath)
            print(new_filepath)

    print("PDB files renamed successfully!")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python pdb_tag.py <pdb_directory>")
    sys.exit(1)

    input_directory = sys.argv[1]
    pdb_tag(input_directory)
