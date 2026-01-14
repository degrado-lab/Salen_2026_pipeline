import sys
import pandas as pd
import os
import shutil
import csv
"""
This script filters the output analysis csv from a lmpnn_rsfr pipeline run and extracts the rosetta model/structure predictions that pass the filters
"""
def filter_dataframe(df, cutoffs):
  """Filters a pandas DataFrame based on given cutoffs for different columns.

  Args:
      df: The pandas DataFrame to be filtered.
      cutoffs: A dictionary where keys are column names and values are tuples of
        (lower_cutoff, upper_cutoff).

  Returns:
      A new pandas DataFrame containing the filtered rows.
  """

  filters = []
  for col, (lower, upper) in cutoffs.items():
    # Handle missing cutoffs and non-numeric columns
    if lower is None and upper is None:
      continue
    elif not pd.api.types.is_numeric_dtype(df[col]):
      raise ValueError(f"Column '{col}' is not numeric. Filtering by cutoffs is not supported.")

    # Create filter condition based on cutoffs
    filters.append((df[col] >= lower) & (df[col] <= upper))

  # Combine filters using logical AND
  combined_filter = filters[0]
  for i in range(1, len(filters)):
    combined_filter = combined_filter & filters[i]

  # Apply combined filter to the dataframe
  return df[combined_filter]

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

def filter_pdbs(df_csv_file, plddt_low, rmsd_high, charge_low, charge_high):
  """
  Sets up inputs for filter_dataframe function, runs filter_dataframe function, outputs filtered CSV

  Args:
      df_csv_file: full run analysis csv file
      plddt_low: minimum plddt cutoff
      rmsd_high: maximum rmsd cutoff 
      charge_low/charge_high: net charge range

  Returns:
      filtered CSV file
  """
  pd_df = pd.read_csv(df_csv_file)
  print('csv read to pandas dataframe')
  cutoffs = {'pLDDT': (plddt_low, 100), 
              'RMSD': (0, rmsd_high),
              'Net charge': (charge_low, charge_high)}  
  filtered_df = filter_dataframe(pd_df.copy(), cutoffs)
  filtered_df_rounded = filtered_df.round(2)

  # Export the rounded DataFrame to CSV
  filtered_df_rounded.to_csv('lmpnn_rsfr_rxs_analysis_filtered.csv', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python filter_pdbs.py csv_file_to_filter plddt_low rmsd_high charge_low charge_high")
        sys.exit(1)

    csv_file = sys.argv[1]
    plddt_low = float(sys.argv[2])
    rmsd_high = float(sys.argv[3])
    charge_low = int(sys.argv[4])
    charge_high = int(sys.argv[5])

    filter_pdbs(csv_file, plddt_low, rmsd_high, charge_low, charge_high)
