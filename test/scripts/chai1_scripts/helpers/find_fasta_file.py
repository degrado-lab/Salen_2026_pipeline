import os
import sys

def find_fasta_file(directory, search_string):
    """
    Searches for files in the given directory with filenames containing the specified string.

    Args:
        directory (str): The directory to search in.
        search_string (str): The string to search for in filenames.

    Returns:
        str: The full path to the first matching file, or None if no match is found.
    """

    for filename in os.listdir(directory):
        if search_string in filename:
            file_path = os.path.join(directory, filename)
            print (file_path) #this print statement is captured by the bash script variable

    return None  # No match found

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(1)

    input_directory = sys.argv[1]
    pdb_tag_id = sys.argv[2]
    find_fasta_file(input_directory, pdb_tag_id)