# Load required libraries
library(sf)
library(dplyr)

# Input and output file paths
centroids_file <- snakemake@input[["merged_file"]]
output_file <- snakemake@output[["pop_centroids"]]

# Load the data
data <- read.csv(centroids_file)

# Check if the required columns exist
if (!all(c("FID", "longitude", "latitude") %in% colnames(data))) {
  stop("The merged file does not contain the required columns: FID, longitude, and latitude.")
}

# Convert data to an sf object, using the correct column names for longitude and latitude
data_sf <- st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326)

# Calculate centroids for each population (FID)
centroids <- data_sf %>%
  group_by(FID) %>%  # Group by population
  summarize(geometry = st_centroid(st_combine(geometry))) %>%
  mutate(longitude = st_coordinates(geometry)[, 1],
         latitude = st_coordinates(geometry)[, 2])

# Write the result to a CSV file
centroids %>%
  st_drop_geometry() %>%  # Drop the geometry column, keeping only lon/lat and FID
  write.csv(output_file, row.names = FALSE)
