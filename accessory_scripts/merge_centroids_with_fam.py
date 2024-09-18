import pandas as pd

# File paths
centroids_file = snakemake.input.ind_centroids
fam_file = snakemake.input.fam
output_file = snakemake.output.merged_file

# Load the individual centroids CSV file
centroids_data = pd.read_csv(centroids_file)

# Load the .fam file (assuming it's space-delimited)
fam_data = pd.read_csv(fam_file, sep='\s+', header=None, names=['FID', 'IndividualID', 'PaternalID', 'MaternalID', 'Sex', 'Phenotype'])

# Merge the two dataframes based on IndividualID
merged_data = centroids_data.merge(fam_data[['IndividualID', 'FID']], on='IndividualID', how='left')

# Reorder the columns to prepend 'FID' (popID)
cols = ['FID'] + [col for col in merged_data.columns if col != 'FID']
merged_data = merged_data[cols]

# Save the merged data to a new CSV file
merged_data.to_csv(output_file, index=False)
