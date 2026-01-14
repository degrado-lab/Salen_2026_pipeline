import os
import sys

def rename_raptor_x_single(directory):
    """Renames PDB files in the given directory, removing the extra extension."""
    os.chdir(directory)  # Change working directory to the specified directory
    for filename in os.listdir():
        if filename.endswith(".RaptorX-Single-ESM1b-ESM1v-ProtTrans.pdb"):
            new_filename = filename[:-41] + ".pdb"  # Remove the extra extension
            os.rename(filename, new_filename)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(1)

    input_directory = sys.argv[1]
    rename_raptor_x_single(input_directory)


