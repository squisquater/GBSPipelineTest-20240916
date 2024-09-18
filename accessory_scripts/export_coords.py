import pandas as pd

# Load the CSV file
df = pd.read_csv(snakemake.input.centroids)

# Extract the 'longitude' and 'latitude' columns
coords_df = df[['longitude', 'latitude']]

# Save the coordinates as a tab-delimited text file without headers
coords_df.to_csv(snakemake.output.coord, sep=' ', index=False, header=False)

print(f"Coordinates saved to {snakemake.output.coord} as a space-delimited file without headers.")
