# mds_plot.R

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(viridis)

# Read the .mds file
mds_data <- read.table(snakemake@input[['mds']], header = TRUE)

# Use the FID column to represent populations
mds_data$Population <- mds_data$FID

# Plotting the MDS
mds_plot <- ggplot(data = mds_data, aes(x = C1, y = C2)) +
  geom_point(aes(color = Population)) +
  labs(title = "MDS Plot",
       x = "C1",
       y = "C2") +
  scale_color_hue()

# Save the plot
ggsave(filename = snakemake@output[['mds_plot']], plot = mds_plot, width = 10, height = 6)
