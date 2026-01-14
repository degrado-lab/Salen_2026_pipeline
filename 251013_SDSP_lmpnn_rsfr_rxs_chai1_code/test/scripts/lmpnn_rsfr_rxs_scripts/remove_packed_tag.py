import os
import re
import sys

def removed_packed_tag(directory):
  """Renames PDB files in a given directory with a specific naming pattern.

  Args:
    directory: The directory containing the PDB files.
  """

  pattern = r"(.+?)_packed_(\d+)_(\d+)\.pdb"

  for filename in os.listdir(directory):
    if filename.endswith('.pdb'):
      match = re.match(pattern, filename)
      if match:
        prefix, numberX, numberY = match.groups()
        new_filename = prefix + "_" + numberX + ".pdb"
        os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))

if __name__ == '__main__':
  if len(sys.argv) != 2: 
    sys.exit(1)
  directory = sys.argv[1]
  removed_packed_tag(directory)