import pandas as pd
import re

# Get the input and output file paths from Snakemake
ped_file = snakemake.input.ped
map_file = snakemake.input.map
output_file = snakemake.output.genhet
locus_output_file = snakemake.output.loci  # Locus names output file

# Function to convert alleles to numeric values
def allele_to_numeric(allele):
    allele_map = {'A': '1', 'C': '2', 'T': '3', 'G': '4', '0': 'NA'}
    return allele_map.get(allele, allele)  # Default to original value if not found

def clean_locus_name(locus):
    """Remove unwanted characters from locus names."""
    # Use regular expressions to remove ':', '+' and '-' characters
    return re.sub(r'[:+/-]', '_', locus)

def plink_to_genhet(ped_file, map_file, output_file, locus_output_file):
    # Read .ped file (individuals and genotype data)
    ped_data = pd.read_csv(ped_file, delim_whitespace=True, header=None)
    
    # Read .map file (SNP information)
    map_data = pd.read_csv(map_file, delim_whitespace=True, header=None, names=["chr", "snp", "cm", "pos"])
    
    # Extract individual IDs (First two columns are Family ID and Individual ID)
    individuals = ped_data.iloc[:, 1]  # The second column contains individual identifiers
    
    # Extract genotype data (the rest of the columns contain genotypes: two columns per locus)
    genotype_data = ped_data.iloc[:, 6:]  # Genotype starts from the 7th column
    
    # Convert genotype data to strings (to ensure consistent replacement)
    genotype_data = genotype_data.astype(str)
    
    # Apply allele to numeric conversion and replace '0' with 'NA'
    genotype_data = genotype_data.applymap(allele_to_numeric)
    
    # Extract and clean locus names from the .map file (2nd column)
    loci_names = map_data['snp'].apply(clean_locus_name).tolist()  # Cleaned locus names
    
    # Save the cleaned locus names to a file
    with open(locus_output_file, 'w') as f:
        for locus in loci_names:
            f.write(f"{locus}\n")
    
    # Create a list to store columns for each locus
    genhet_columns = []
    
    # Loop through the genotype columns two at a time (each locus has two columns)
    for i in range(0, genotype_data.shape[1], 2):
        locus_alleles = genotype_data.iloc[:, i:i+2]
        locus_alleles.columns = [f'Locus_{i//2+1}_A1', f'Locus_{i//2+1}_A2']  # Temporary column names
        genhet_columns.append(locus_alleles)
    
    # Concatenate all locus columns horizontally (axis=1)
    genotype_concat = pd.concat(genhet_columns, axis=1)
    
    # Create the final GenHet input DataFrame without headers
    genhet_data = pd.DataFrame()
    genhet_data['IndividualID'] = individuals
    genhet_data = pd.concat([genhet_data, genotype_concat], axis=1)
    
    # Write to output file without header (tab-separated format)
    genhet_data.to_csv(output_file, sep='\t', index=False, header=False)

# Call the function
plink_to_genhet(ped_file, map_file, output_file, locus_output_file)
