import os
import sys 

def mpnn_out_to_single_file_fasta(pdb_directory, output_directory):
    # Output fasta file name
    output_file = f"{output_directory}/single_file.fasta"

    # List to store sequences
    sequences = []

    # Loop through each file in the directory
    for filename in os.listdir(pdb_directory):
        if filename.endswith(".fa"):
            print (filename)
            # Open the file for reading
            with open(os.path.join(pdb_directory, filename), "r") as f:
                # Read in the contents of the file
                contents = f.read()

            # Split the contents of the file into individual sequences
            seq_info = contents.split(">")[1:]

            # Flag to skip the first sequence
            skip_first = True

            # Loop through each sequence and modify the header
            for seq in seq_info:
                # Skip the first sequence
                if skip_first:
                    skip_first = False
                    continue

                # Split the sequence information from the sequence itself
                header, sequence = seq.split("\n", 1)

                # Extract the sample number from the header
                sample = header.split(",")[1].split("=")[1]

                # Create a new header using the file name and sample number
                new_header = filename.split(".")[0] + "_" + sample

                # Append the modified sequence to the list
                sequences.append(">" + new_header + "\n" + sequence + "\n")

    # Open the output file for writing
    with open(output_file, "w") as f:
        # Write all the sequences to the output file
        for seq in sequences:
            f.write(seq)

    print("All sequences have been written to:", output_file)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(1)

    input_directory = sys.argv[1]
    output_directory = sys.argv[2]
    mpnn_out_to_single_file_fasta(input_directory, output_directory)

