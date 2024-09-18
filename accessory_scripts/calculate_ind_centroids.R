# Load required libraries
library(sf)
library(dplyr)

# Input and output file paths from Snakemake
centroids_file <- snakemake@input[["ind_coords"]]
output_file <- snakemake@output[["ind_centroids"]]

# Debugging: Print the input and output file paths
print(paste("Input file:", centroids_file))
print(paste("Output file:", output_file))

# Load the data
data <- read.csv(centroids_file)

# Convert any non-numeric values in the 'Long' and 'Lat' columns
data <- data %>%
  mutate(Long = as.numeric(Long), Lat = as.numeric(Lat))

# Log rows with issues in Long or Lat
problematic_rows <- data %>%
  filter(is.na(Long) | is.na(Lat))

if (nrow(problematic_rows) > 0) {
  write.csv(problematic_rows, file="problematic_coords.csv", row.names = FALSE)
  message("Warning: Found rows with missing or invalid coordinates. See 'problematic_coords.csv' for details.")
}

# Convert data to an sf object, assuming longitude and latitude columns are named 'Long' and 'Lat'
data_sf <- st_as_sf(data, coords = c("Long", "Lat"), crs = 4326)

# Calculate centroids for replicate samples
centroids <- data_sf %>%
  group_by(IndividualID) %>%  # Assuming "IndividualID" is the column name for sample IDs
  summarize(geometry = st_centroid(st_combine(geometry))) %>%
  mutate(longitude = st_coordinates(geometry)[, 1],
         latitude = st_coordinates(geometry)[, 2])

# Write the result to a CSV file
centroids %>%
  st_drop_geometry() %>%  # Drop the geometry column, keeping only lon/lat and sample ID
  write.csv(output_file, row.names = FALSE)
