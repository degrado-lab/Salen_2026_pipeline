import os
import re
import sys

def mpnn_out_to_raptor_x_input(pdb_directory, output_directory): 
    # New output folder name
    output_folder = f"{output_directory}/fastas"

    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Loop through each file in the directory
    for filename in os.listdir(pdb_directory):
        if filename.endswith(".fa"):
            # Open the file for reading
            with open(os.path.join(pdb_directory, filename), "r") as f:
                # Read in the contents of the file
                contents = f.read()

            # Split the contents of the file into individual sequences
            seq_info = contents.split(">")[1:]

            # Skip the first sequence for each file
            if seq_info:
                seq_info = seq_info[1:]  # Discard the first sequence

            # Loop through each remaining sequence and create a separate file
            for seq in seq_info:
                # Split the sequence information from the sequence itself
                header, sequence = seq.split("\n", 1)

                # Extract ID number (handle potential missing comma)
                id_number = ""
                if ", id=" in header:
                    id_number = header.split(", id=")[1].split(", ")[0]

                # Create new header using the file name and ID number (or default)
                new_header = f"{filename.split('.')[0]}_{id_number}" if id_number else filename.split('.')[0]

                # Create the output file path
                output_filename = os.path.join(output_folder, f"{new_header}.fasta")

                # Write the modified sequence to the output file
                with open(output_filename, "w") as f:
                    f.write(f">{new_header}\n{sequence}\n")

    print("All sequences (except the first one in each file) have been written to individual files in:", output_folder)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(1)

    input_directory = sys.argv[1]
    output_directory = sys.argv[2]
    mpnn_out_to_raptor_x_input(input_directory, output_directory)
