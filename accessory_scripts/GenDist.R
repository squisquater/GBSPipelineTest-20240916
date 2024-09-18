# Load necessary libraries
library(adegenet)
library(vcfR)
library(hierfstat)
library(StAMPP)
library(ggplot2)
library(viridis)

# Pull inputs directly from Snakemake
landgen_dir <- snakemake@params[["landgen_dir"]]
vcf_file <- snakemake@input[["vcf"]]
fam_file <- snakemake@input[["fam"]]

# Set working directory to the LandGen directory
setwd(landgen_dir)

# Read the population information file from the FAM file
fam_data <- read.table(fam_file, header = FALSE)
pop_data <- data.frame(
  Population = fam_data$V1,
  IndividualID = paste(fam_data$V1, fam_data$V2, sep = "_")
)

# Read the VCF file
vcf <- read.vcfR(vcf_file)

# Convert VCF to genlight object
genlight_obj <- vcfR2genlight(vcf)

# Debugging step: check the sample IDs in the genlight object
print(genlight_obj@ind.names)

# Ensure the sample IDs in pop_data match those in genlight_obj
if (!all(pop_data$IndividualID %in% genlight_obj@ind.names)) {
  stop("Some sample IDs in FAM file do not match those in the VCF file.")
} else {
  message("All sample IDs in FAM file match those in the VCF file.")
}

# Assign population information to genlight object
pop(genlight_obj) <- pop_data$Population[match(genlight_obj@ind.names, pop_data$IndividualID)]

# Verify the populations have been correctly assigned
print(pop(genlight_obj))

# Convert genlight object to StAMPP format
stampp_data <- stamppConvert(genlight_obj, "genlight")

# Calculate pairwise Fst
stampp_fst <- stamppFst(stampp_data, nboots = 100, percent = 95, nclusters = 1)

# Extract the Fst matrix
fst_matrix <- stampp_fst$Fsts

# Fill in the entire Fst matrix by mirroring the lower triangle to the upper triangle
fst_matrix[upper.tri(fst_matrix)] <- t(fst_matrix)[upper.tri(fst_matrix)]

# Replace NA values with 0
fst_matrix[is.na(fst_matrix)] <- 0

# Print the Fst matrix
print(fst_matrix)

# Save the Fst matrix to a CSV file
write.csv(fst_matrix, file.path(landgen_dir, "pairwise_fst_matrix.csv"), row.names = TRUE)

# Calculate Nei's genetic distance (D)
neis_d_matrix <- stamppNeisD(stampp_data)

# Convert the result to a matrix
neis_d_matrix <- as.matrix(neis_d_matrix)

# Assign population IDs as row and column names
colnames(neis_d_matrix) <- rownames(neis_d_matrix)

# Print Nei's D matrix
print(neis_d_matrix)

# Save Nei's D matrix to a CSV file
write.csv(neis_d_matrix, file.path(landgen_dir, "neis_d_matrix.csv"), row.names = TRUE)

#################################
##### Plotting the Matrices #####
#################################

# Convert the Fst matrix to a format suitable for ggplot
fst_df <- as.data.frame(as.table(fst_matrix))

# Plot the pairwise Fst matrix as a heatmap
fst_heatmap <- ggplot(fst_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_viridis(option = "magma", name = "Pairwise Fst", direction = -1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Population", y = "Population", title = "Pairwise Fst Heatmap")

# Save the Fst heatmap
ggsave(file.path(landgen_dir, "pairwise_fst_heatmap.png"), fst_heatmap, width = 10, height = 8)

# Convert the Nei's D matrix to a format suitable for ggplot
neis_d_df <- as.data.frame(as.table(neis_d_matrix))

# Plot the pairwise Nei's D matrix as a heatmap
neis_d_heatmap <- ggplot(neis_d_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_viridis(option = "magma", name = "Pairwise Nei's D", direction = -1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Population", y = "Population", title = "Pairwise Nei's D Heatmap")

# Save the Nei's D heatmap
ggsave(file.path(landgen_dir, "pairwise_neis_d_heatmap.png"), neis_d_heatmap, width = 10, height = 8)
