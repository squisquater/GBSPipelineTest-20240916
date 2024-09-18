#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)

# Input and output file paths
genhet_file <- snakemake@input$genhet
depth_file <- snakemake@input$depth
fam_file <- snakemake@input$fam
output_file <- snakemake@output$merged

# Load the genetic heterozygosity metrics file
genhet_data <- read.table(genhet_file, header = TRUE, sep = "\t", na.strings = "NA")

# Remove ".merged" suffix from sample IDs in the genhet file
genhet_data$sampleid <- gsub("\\.merged$", "", genhet_data$sampleid)

# Load the depth file
depth_data <- read.table(depth_file, header = TRUE, sep = "\t")

# Remove duplicates from depth file
depth_data <- depth_data %>% distinct(BamFile, .keep_all = TRUE)

# Load the .fam file and extract the Region.ID
fam_data <- read.table(fam_file, header = FALSE, sep = "")
colnames(fam_data) <- c("FamilyID", "sampleid", "paternalID", "maternalID", "sex", "Region.ID")

# Remove ".merged" suffix from sample IDs in the fam file
fam_data$sampleid <- gsub("\\.merged$", "", fam_data$sampleid)

# Select FamilyID as Region.ID from fam_data
fam_info <- fam_data %>%
  select(sampleid, FamilyID) %>%
  rename(Region.ID = FamilyID)

# Merge the depth, genetic metrics, and population information
merged_data <- genhet_data %>%
  inner_join(depth_data, by = c("sampleid" = "BamFile")) %>%
  inner_join(fam_info, by = "sampleid")

# Save the merged file
write.table(merged_data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Merged data with population info saved to:", output_file, "\n")
