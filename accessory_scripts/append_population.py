import pandas as pd
import argparse

def prepend_population_info(vcf_file, popmap_file, output_vcf):
    # Load the popmap
    popmap = pd.read_csv(popmap_file, sep="\t", header=None, names=["SampleID", "Population"])
    
    # Read the VCF file
    with open(vcf_file, 'r') as vcf_in:
        vcf_lines = vcf_in.readlines()
    
    with open(output_vcf, 'w') as vcf_out:
        for line in vcf_lines:
            if line.startswith("#CHROM"):
                # Modify the header line to prepend population info to each sample ID
                headers = line.strip().split("\t")
                sample_columns = headers[9:]  # Sample columns start from the 10th column
                new_columns = []
                for sample in sample_columns:
                    pop = popmap.loc[popmap['SampleID'] == sample, 'Population'].values
                    if len(pop) == 0:
                        pop = "Unknown"  # Handle cases where the sample is not found in the popmap
                    else:
                        pop = pop[0]
                    # Prepend the population to the SampleID
                    new_columns.append(f"{pop}_{sample}")
                vcf_out.write("\t".join(headers[:9] + new_columns) + "\n")
            else:
                vcf_out.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepend population information to VCF file.")
    parser.add_argument("--input", required=True, help="Input VCF file")
    parser.add_argument("--output", required=True, help="Output VCF file with population info prepended")
    parser.add_argument("--popmap", required=True, help="Popmap file with SampleID and Population")
    
    args = parser.parse_args()
    
    prepend_population_info(args.input, args.popmap, args.output)
