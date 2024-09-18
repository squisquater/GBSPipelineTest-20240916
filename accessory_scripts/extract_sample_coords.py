import pandas as pd

# File paths
fam_file = snakemake.input.fam
database_file = snakemake.input.database
output_file = snakemake.output.ind_coords

# Load the .fam file (assuming it's space-delimited and individual IDs are in the second column)
fam_data = pd.read_csv(fam_file, sep='\s+', header=None, usecols=[1], names=['IndividualID'])

# Create a new column without the ".merged" suffix for matching
fam_data['IndividualID_no_suffix'] = fam_data['IndividualID'].str.replace('.merged$', '', regex=True)

# Load the master database file
database = pd.read_csv(database_file, sep='\t')  # Assuming it's tab-delimited

# Select the relevant columns from the database
database_subset = database[['Individual.ID', 'Lat', 'Long']]

# Merge based on matching Individual IDs (without .merged suffix)
merged_data = fam_data.merge(database_subset, left_on='IndividualID_no_suffix', right_on='Individual.ID')

# Save the output, keeping the original IndividualID (with .merged)
merged_data[['IndividualID', 'Lat', 'Long']].to_csv(output_file, index=False)
