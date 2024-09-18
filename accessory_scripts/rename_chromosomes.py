import pandas as pd
import sys

# Command line arguments
bim_file_path = sys.argv[1]
output_bim_file_path = sys.argv[2]
renaming_key_path = sys.argv[3]

# Load .bim file
bim_data = pd.read_csv(bim_file_path, sep='\t', header=None, names=["chrom", "snp", "gen_dist", "bp_position", "allele1", "allele2"])

# Generate renaming scheme and apply it
unique_chroms = bim_data["chrom"].unique()
renaming_scheme = {chrom: str(i+1) for i, chrom in enumerate(sorted(unique_chroms))}
bim_data["chrom"] = bim_data["chrom"].map(renaming_scheme)

# Save modified .bim file
bim_data.to_csv(output_bim_file_path, sep='\t', index=False, header=False)

# Save renaming key
with open(renaming_key_path, 'w') as f:
    for key, value in renaming_scheme.items():
        f.write(f"{key}\t{value}\n")
