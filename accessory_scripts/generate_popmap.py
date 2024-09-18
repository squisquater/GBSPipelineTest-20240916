import pandas as pd

# Snakemake provides input and output as variables
depth_file = str(snakemake.input.depth_filtered)  # Depth-filtered file from Snakemake
database_file = str(snakemake.input.database)  # Database file from Snakemake
output_file = str(snakemake.output.popmap)  # Output popmap file from Snakemake

# Read the depth-filtered file (samples with the ".merged" suffix)
depth_df = pd.read_csv(depth_file, header=None, names=["SampleID"])
depth_df['SampleID'] = depth_df['SampleID'].str.replace('.merged', '')  # Remove ".merged" to match database

# Read the database file
db_df = pd.read_csv(database_file, sep="\t")

# Merge the two based on matching Individual.ID from the database and SampleID from depth file
merged_df = pd.merge(depth_df, db_df[['Individual.ID', 'Region.ID']], left_on="SampleID", right_on="Individual.ID", how="left")

# Add back the ".merged" suffix to the SampleID
merged_df['SampleID'] = merged_df['SampleID'] + ".merged"

# Remove duplicates
merged_df = merged_df.drop_duplicates(subset=['SampleID', 'Region.ID'])

# Sort by 'Region.ID' (Population) first, then by 'SampleID'
merged_df = merged_df.sort_values(by=['Region.ID', 'SampleID'])

# Save the result without headers and index
merged_df[['SampleID', 'Region.ID']].to_csv(output_file, sep="\t", index=False, header=False)
