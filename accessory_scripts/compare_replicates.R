library(dplyr)

# Load in the output file of your IBS analysis
main_data <- read.table(snakemake@input[['ibs_file']], header = TRUE)

# Load in your sampleID key (this should have numbers 1-X in col1 and associated sample ID in col2)
id_mapping <- read.table(snakemake@input[['replicate_key']], header = FALSE)
colnames(id_mapping) <- c("FID", "replacement_ID")

# Load the library depth summary
depth_data <- read.table(snakemake@input[['depth_file']], header = TRUE)  # Assuming 'depth_file' is your path to the library_depth_summary.txt

# Replace FID1 and FID2 with meaningful IDs
main_data <- main_data %>%
  left_join(id_mapping, by = c("FID1" = "FID")) %>%
  mutate(FID1 = ifelse(!is.na(replacement_ID), replacement_ID, FID1)) %>%
  select(-replacement_ID)

main_data <- main_data %>%
  left_join(id_mapping, by = c("FID2" = "FID")) %>%
  mutate(FID2 = ifelse(!is.na(replacement_ID), replacement_ID, FID2)) %>%
  select(-replacement_ID)

# Initialize the results dataframe
results <- data.frame(
  Individual = character(),
  IndividualReadDepth = numeric(),
  MaxPartner = character(),
  MaxPartnerReadDepth = numeric(),
  MaxPI_HAT = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each individual
for (individual in unique(c(main_data$FID1, main_data$FID2))) {
  individual_data <- main_data %>% 
    filter(FID1 == individual | FID2 == individual) %>%
    filter(!is.nan(PI_HAT)) %>%
    filter(PPC > 0.99)
  
  if(nrow(individual_data) > 0) {
    max_row <- individual_data[which.max(individual_data$PI_HAT),]
    max_partner <- ifelse(max_row$FID1 == individual, max_row$FID2, max_row$FID1)
    
    # Get read depth for individual and max partner
individual_index <- match(individual, depth_data$BamFile)
max_partner_index <- match(max_partner, depth_data$BamFile)

individual_depth <- if(!is.na(individual_index)) depth_data$GenomeWideAvgDepth[individual_index] else NA
max_partner_depth <- if(!is.na(max_partner_index)) depth_data$GenomeWideAvgDepth[max_partner_index] else NA

# Append the result
results <- rbind(results, data.frame(
  Individual = individual,
  IndividualReadDepth = individual_depth,
  MaxPartner = max_partner,
  MaxPartnerReadDepth = max_partner_depth,
  MaxPI_HAT = max_row$PI_HAT
))
  } else {
    # Handle cases with no valid comparisons
    results <- rbind(results, data.frame(
      Individual = individual,
      IndividualReadDepth = NA,
      MaxPartner = NA,
      MaxPartnerReadDepth = NA,
      MaxPI_HAT = NA
    ))
  }
}

# Write the results to a .csv file
write.csv(results, snakemake@output[['pi_hat_replicates']], row.names = FALSE)
