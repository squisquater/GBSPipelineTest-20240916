library(adegenet)
library(vcfR)
library(hierfstat)
library(StAMPP)
library(ggplot2)
library(viridis)

# Set working directory to the LandGen directory defined in Snakemake
landgen_dir <- snakemake@params[['landgen_dir']]
setwd(landgen_dir)

# Read the population information file from the fam file
fam_file <- read.table(snakemake@input[['fam']], header = FALSE)
pop_data <- data.frame(
  Population = fam_file$V1,
  IndividualID = paste(fam_file$V1, fam_file$V2, sep = "_")
)

# Read the VCF file
vcf <- read.vcfR(snakemake@input[['vcf']])

# Convert VCF to a genind object
genind_obj <- vcfR2genind(vcf)

# Check if the genind object is created properly
if (!inherits(genind_obj, "genind")) {
  stop("The object is not a valid genind object.")
}

# Access individual names from the row names of the tab matrix
ind_names <- rownames(genind_obj@tab)
print(ind_names)

# Ensure the sample IDs in pop_data match those in the genind object
if (!all(pop_data$IndividualID %in% ind_names)) {
  stop("Some sample IDs in populations_info.txt do not match those in the VCF file.")
} else {
  message("All sample IDs in populations_info.txt match those in the VCF file.")
}

# Assign population information to genind object
pop(genind_obj) <- pop_data$Population[match(ind_names, pop_data$IndividualID)]

# Verify the populations have been correctly assigned
print(pop(genind_obj))

# Convert genind object to a genpop object
genpop_obj <- genind2genpop(genind_obj)

# Ensure the subpopulation IDs in geo_centroids match those in genpop_obj
genpop_populations <- row.names(tab(genpop_obj))

# Read the geographic information file for subpopulation centroids (Ensure columns FID, longitude, latitude)
geo_centroids <- read.csv(snakemake@input[['pop_centroids']])

# Print the column names of geo_centroids to confirm correct structure
print(colnames(geo_centroids))

# Check if all required columns exist
if (!all(c("FID", "longitude", "latitude") %in% colnames(geo_centroids))) {
  stop("The population centroids file does not contain the required columns: FID, longitude, and latitude.")
}

# Ensure the subpopulation IDs in geo_centroids match those in genpop_obj
if (!all(genpop_populations %in% geo_centroids$FID)) {
  stop("Some population IDs in population_centroids.csv do not match those in the genpop object.")
} else {
  message("All population IDs in population_centroids.csv match those in the genpop object.")
}

# Assign latitude and longitude information to a matrix using the correct column names
latlong <- geo_centroids[match(genpop_populations, geo_centroids$FID), c("longitude", "latitude")]
rownames(latlong) <- genpop_populations

# Print the latlong matrix to ensure correct extraction
print(latlong)

# Calculate pairwise genetic distances and pairwise geographic distances
Dgen <- dist.genpop(genpop_obj)
Dgeo <- dist(latlong)

# Check for missing values in the matrices
if (any(is.na(Dgeo))) {
  stop("Geographic distance matrix contains missing values.")
}
if (any(is.na(Dgen))) {
  stop("Genetic distance matrix contains missing values.")
}

# Full Mantel test (with all populations)
ibd_result_full <- mantel.randtest(Dgen, Dgeo)
print(ibd_result_full)

# Extract relevant information from the full Mantel test result
mantel_result_data_full <- data.frame(
  obs = ibd_result_full$obs,      
  pvalue = ibd_result_full$pvalue,
  alternative = ibd_result_full$alter,
  replicates = ibd_result_full$rep
)

# Save the full Mantel test result to a text file
write.table(mantel_result_data_full, file = snakemake@output[['mantel_result']], quote = FALSE, row.names = FALSE, col.names = TRUE)

# Subset populations and run Mantel test excluding listed populations (if any)
exclude_pops <- snakemake@params[['exclude_pops']]

# Subset the populations to exclude from the matrices
keep_pops <- setdiff(genpop_populations, exclude_pops)

Dgen_exclude <- as.dist(as.matrix(Dgen)[keep_pops, keep_pops])
Dgeo_exclude <- as.dist(as.matrix(Dgeo)[keep_pops, keep_pops])

# Run Mantel test with excluded populations (or full dataset if exclude_pops is empty)
ibd_result_exclude <- mantel.randtest(Dgen_exclude, Dgeo_exclude)

# Extract relevant information from the exclude Mantel test result
mantel_result_data_exclude <- data.frame(
  obs = ibd_result_exclude$obs,
  pvalue = ibd_result_exclude$pvalue,
  alternative = ibd_result_exclude$alter,
  replicates = ibd_result_exclude$rep
)

# Save the exclude Mantel test result to a separate text file
write.table(mantel_result_data_exclude, file = snakemake@output[['mantel_result_exclude']], quote = FALSE, row.names = FALSE, col.names = TRUE)

# Export geographic and genetic distance matrices (as before)
write.csv(as.matrix(Dgeo), file = snakemake@output[['geo_dist_matrix']], row.names = TRUE)
write.csv(as.matrix(Dgen), file = snakemake@output[['gen_dist_matrix']], row.names = TRUE)

# Convert dist objects to matrices for full Mantel test
geo_matrix_full <- as.matrix(Dgeo)
gen_matrix_full <- as.matrix(Dgen)

# Extract the lower triangle of the matrices as vectors for full Mantel test
geo_vector_full <- geo_matrix_full[lower.tri(geo_matrix_full)]
gen_vector_full <- gen_matrix_full[lower.tri(gen_matrix_full)]

# Fit a linear model for full Mantel test
model_full <- lm(gen_vector_full ~ geo_vector_full)

# Plot geographic vs. genetic distances with the linear model for full Mantel test
plot_data_full <- data.frame(GeographicDistance = geo_vector_full, GeneticDistance = gen_vector_full)

# Debugging: print a summary of the linear model
print(summary(model_full))

# Save the plot for the full Mantel test to a PNG file
png(snakemake@output[['ibd_plot']], width = 800, height = 600)
ggplot(plot_data_full, aes(x = GeographicDistance, y = GeneticDistance)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red") +
  labs(title = paste("Mantel Test: r =", round(ibd_result_full$obs, 3), "p =", round(ibd_result_full$pvalue, 3)),
       x = "Geographic distance", y = "Genetic distance") +
  annotate("text", x = Inf, y = Inf, label = paste("r =", round(cor(geo_vector_full, gen_vector_full), 3)), hjust = 1.1, vjust = 1.1)
dev.off()

# Always create a plot for ibd_result_exclude
geo_matrix_exclude <- as.matrix(Dgeo_exclude)
gen_matrix_exclude <- as.matrix(Dgen_exclude)

# Extract the lower triangle of the matrices as vectors for the exclude test
geo_vector_exclude <- geo_matrix_exclude[lower.tri(geo_matrix_exclude)]
gen_vector_exclude <- gen_matrix_exclude[lower.tri(gen_matrix_exclude)]

# Fit a linear model for the exclude test
model_exclude <- lm(gen_vector_exclude ~ geo_vector_exclude)

# Plot geographic vs. genetic distances with the linear model for the exclude test
plot_data_exclude <- data.frame(GeographicDistance = geo_vector_exclude, GeneticDistance = gen_vector_exclude)

# Save the plot for the Mantel test with excluded populations
png(snakemake@output[['ibd_plot_exclude']], width = 800, height = 600)
ggplot(plot_data_exclude, aes(x = GeographicDistance, y = GeneticDistance)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue") +
  labs(title = paste("Mantel Test (Exclude Populations): r =", round(ibd_result_exclude$obs, 3), "p =", round(ibd_result_exclude$pvalue, 3)),
       x = "Geographic distance", y = "Genetic distance") +
  annotate("text", x = Inf, y = Inf, label = paste("r =", round(cor(geo_vector_exclude, gen_vector_exclude), 3)), hjust = 1.1, vjust = 1.1)
dev.off()