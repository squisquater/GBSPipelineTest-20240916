# Load necessary libraries
library(dplyr)

# Read arguments from the snakemake object
admixture_dir <- snakemake@params$admixture_dir
fam_file <- snakemake@params$fam_file
q_files <- snakemake@input$q_files  # Directly access the q_files passed from Snakemake

# Define the output file path within the admixture_dir
output_file <- file.path(admixture_dir, "merged_q_values_summary.txt")

# Function to generate column names based on K
generate_colnames <- function(K) {
  paste0("K", K, "_Q", seq_len(K))
}

# Initialize list to store data frames
q_data_list <- list()

# Loop over the .Q files and read them
for (q_file in q_files) {
  # Extract K value from the file name
  K <- as.numeric(gsub(".*_renamedChr\\.(\\d+)\\.Q$", "\\1", q_file))
  
  # Read the Q file
  q_data <- read.table(q_file, header = FALSE)
  cat("Read Q file:", q_file, "for K =", K, "with", ncol(q_data), "columns\n")
  
  # Generate column names based on K
  colnames(q_data) <- generate_colnames(K)
  
  # Store the data frame in the list
  q_data_list[[K]] <- q_data
}

# Merge all Q data frames by columns
merged_q_data <- bind_cols(q_data_list)

# Read the .fam file, assuming it has multiple columns
fam_data <- read.table(fam_file, header = FALSE, colClasses = "character")

# Check if fam_data has at least two columns
if (ncol(fam_data) < 2) {
  stop("Error: The .fam file does not have at least two columns as expected.")
}

# Extract the second column (typically the IID) as the individual ID
fam_data <- fam_data[, 2, drop = FALSE]
colnames(fam_data) <- "Individual_ID"

# Join the fam data with the merged Q data
final_data <- cbind(fam_data, merged_q_data)

# Write the output to a file
write.table(final_data, file = output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

cat("Q-values summarized and saved to:", output_file, "\n")
