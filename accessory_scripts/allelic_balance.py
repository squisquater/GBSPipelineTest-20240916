#! /usr/bin/env python

# Author: Joana
# Modified by: SophiePQ on 8/26/2024
# Specifically modified to:
# 1. define sites as monomorphic if AF=0
# 2. use updated binomial test function from scipy.stats
# 3. print the number of sites set to missing due to failing the allelic balance test
# 4. generate a log file summarizing the number of SNPs set to missing per individual

# Usage: python ./allelicBalance_modified.py --input populations.snps.vcf.gz --output AB_output.vcf --log AB_output.log --exclude --pvalue 0.01 --ratio 0.2

from sys import *
import os, time, argparse, re, scipy, numpy
from collections import defaultdict
from scipy.stats import binomtest
import gzip

parser = argparse.ArgumentParser(description='Filter out genotypes with strong allelic disbalance in the given VCF file')

parser.add_argument('-i', '--input', dest='i', help="input file in vcf format [required]", required=True)
parser.add_argument('-o', '--output', dest='o', help="output file [required]", required=True)
parser.add_argument('-l', '--log', dest='log', help="log file [required]", required=True)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-hom', '--homozygote', action='store_true', help="set failing genotypes as homozygous", default=False)
group.add_argument('-excl', '--exclude', action='store_true', help="set failing genotypes as missing", default=False) 
group.add_argument('-two', '--twoSteps', action='store_true', help="set failing genotypes as missing if pvalue<p and as homozygous if pvalue<p2", default=False)
parser.add_argument('-p', '--pvalue', type=float, dest='p', help="p-value threshold for binomial test [default: 0.01]", required=False, default=0.01)
parser.add_argument('-p2', '--pHomoz', type=float, dest='p2', help="second p-value threshold for binomial test if -two is specified, threshold to make failing genotypes homozygote [default=0.005]", default=0.005)
parser.add_argument('-r', '--ratio', type=float, dest='r', help="hard cutoff for allelic ratio [default=0.2]", default=0.2)

args = parser.parse_args()

threshold = args.p
threshold2 = args.p2
ratio = args.r

# Initialize a dictionary to count missing SNPs per individual
missing_counts = defaultdict(int)

# Check if the input file is gzipped
if args.i.endswith(".gz"):
    inputF = gzip.open(args.i, 'rt')
    print("Assuming the input file is gzipped")
else:
    inputF = open(args.i, 'r')
    print("Assuming the input file is not gzipped")

# Prepare the output file (if filename ends with gz, create gzipped file)
if args.o.endswith(".gz"):
    outputF = gzip.open(args.o, 'wt')
    print("Generating a gzipped output file")
else:
    outputF = open(args.o, 'w')
    print("Generating a regular output file")

# Store individual identifiers from the header line
individuals = []

# First, find and store the individual identifiers
line_count = 0
for Line in inputF:
    line_count += 1
    if Line.startswith("#CHROM"):
        columns = Line.strip().split("\t")
        individuals = columns[9:]  # store individual IDs (columns after the 9th column)
        break

# Reset the input file reader after reading the header
inputF.seek(0)

# Process the VCF file
line_count = 0
for Line in inputF:
    line_count += 1
    if line_count <= 10:  # Print the first 10 lines to confirm reading
        print(f"Line {line_count}: {Line.strip()}")

    # DATA SECTION: clause checks if the header section is over
    if re.match('^#', Line) is None:
        print(f"Processing line {line_count}...")

        # Get the columns of that line
        columns = Line.strip("\n").split("\t")
        print(f"Columns: {columns[:9]}")  # Print the first 9 columns for reference
        
        # Get the format field
        format = columns[8].split(":")
        print(f"Format field: {format}")

        # Check if the required AD field is found (reads for each allele)
        if "AD" in format:
            # Get the AD column number
            AD = format.index("AD")
        else:
            outputF.write(Line)
            print(f"Warning: no AD field found at {columns[0]} {columns[1]}")
            continue

        # Use the AF value in INFO to determine if the site is polymorphic
        info_field = columns[7]
        af_value = None
        if "AF=" in info_field:
            af_value = float(re.search(r"AF=([0-9\.]+)", info_field).group(1))

        # Process only if AF is greater than 0 (indicating polymorphism)
        if af_value is not None and af_value > 0:
            print(f"Processing polymorphic SNP at {columns[0]} {columns[1]} with AF={af_value}")

            # Add the info to the site
            result = columns[0:9]

            # Get the genotypes
            genotypecolumns = range(9, len(columns))

            # Check each individual if it is a heterozygote
            for ind in genotypecolumns:
                genotype = columns[ind]
                genotype = genotype.split(":")

                if "/" in genotype[0]:
                    alleles = genotype[0].split("/") 
                elif "|" in genotype[0]:
                    alleles = genotype[0].split("|") 
                else:
                    result.append(":".join(genotype))
                    continue

                # If the genotype is heterozygous check the allelic balance
                if alleles[0] != alleles[1]:
                    reads = genotype[AD].split(",")

                    # Calculate the probability for the observed allele distribution with a binomial test
                    pval = binomtest(int(reads[0]), n=int(reads[0]) + int(reads[1]), p=0.5).pvalue
                    print(f"p-value for genotype {genotype[0]}: {pval}, reads: {reads}")

                    # if one of the alleles has no reads (weirdly this happens)
                    if int(reads[0]) == 0 or int(reads[1]) == 0:
                        result.append("./.")
                        missing_counts[ind] += 1
                        print("Genotype set to missing due to zero reads.")

                    # If significant allelic disbalance
                    elif pval < threshold:
                        print(f"Significant imbalance detected with p-value {pval}")

                        # If -excl is specified, set failing genotypes as missing
                        if args.exclude:
                            result.append("./.")
                            missing_counts[ind] += 1
                            print(f"Genotype set to missing due to significant imbalance.")

                        # if -hom is specified, set failing genotypes as homozygous
                        elif args.homozygote:
                        
                            # Replace the genotype by the more common allele homozygote
                            if int(reads[0]) > int(reads[1]):
                                genotype[0] = alleles[0] + "/" + alleles[0]
                            else:
                                genotype[0] = alleles[1] + "/" + alleles[1]

                            result.append(":".join(genotype))
                            print(f"Genotype set to homozygous: {genotype[0]}")

                        # if -two is specified, set genotypes failing threshold as missing, and those failing threshold2 as homozygous
                        elif args.twoSteps:
                            if pval > threshold2:
                                result.append("./.")
                                missing_counts[ind] += 1
                                print("Genotype set to missing (two-step).")
                            else:
                                if int(reads[0]) > int(reads[1]):
                                    genotype[0] = alleles[0] + "/" + alleles[0]
                                else:
                                    genotype[0] = alleles[1] + "/" + alleles[1]
                                result.append(":".join(genotype))
                                print(f"Genotype set to homozygous (two-step): {genotype[0]}")
                        else: 
                            exit("Either -hom or -excl or -two is required, please specify how failing genotypes should be handled.")

                    # If the binomial test is not significant but the allelic ratio is very small (at low depth possible) 
                    elif float(reads[0]) / float(reads[1]) < ratio or float(reads[1]) / float(reads[0]) < ratio:
                        result.append("./.")
                        missing_counts[ind] += 1
                        print("Genotype set to missing due to allelic ratio.")

                    # If the genotype passes the test
                    else:
                        result.append(":".join(genotype))
                        print("Genotype passed the test.")

                # If the genotype is homozygous, just append as is
                else:
                    result.append(":".join(genotype))
                    print(f"Genotype is homozygous: {genotype[0]}")

            outputF.write('\t'.join(result) + "\n")

        # If it is not a polymorphic site, just write the line out
        else:
            outputF.write(Line)
            print(f"Monomorphic site at {columns[0]} {columns[1]} (AF={af_value}), writing line as is.")

    # If it is a header line, just write it out
    else:
        outputF.write(Line)
        print("Writing header line.")

# Write the log file summarizing the number of SNPs set to missing per individual
with open(args.log, 'w') as logF:
    logF.write("Individual\tMissing SNPs\n")
    for idx, ind in enumerate(sorted(missing_counts)):
        logF.write(f"{individuals[ind-9]}\t{missing_counts[ind]}\n")  # Adjust index to match individuals list
    print(f"Log file written to {args.log}")

print("\nAll done")

inputF.close()
outputF.close()