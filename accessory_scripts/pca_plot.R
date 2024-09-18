# pca_plot.R

# ############## PCA ################

# Load the ggplot2, dplyr, and gridExtra libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(viridis)


# Read the .eigenval and .eigenvec files
eigenvals <- read.table(snakemake@input[['eigenval']])
eigenvecs <- read.table(snakemake@input[['eigenvec']])

populations <- eigenvecs[, 1]

# Plotting the PCA
pca_plot <- ggplot(data = eigenvecs, aes(x = V3, y = V4)) +
  geom_point(aes(color = populations)) +
  labs(title = "PCA Plot",
       x = "PC1",
       y = "PC2") +
  scale_color_hue()
  #scale_color_manual(values = c("red", "blue", "green")) # You can set your own colors here if you prefer

# Calculate the sum of all eigenvalues
total_eigenvalue <- sum(eigenvals$V1)

# Compute the percentage for each eigenvalue
eigenvals$percentage <- (eigenvals$V1 / total_eigenvalue) * 100

# Subset to only keep the first 4 eigenvalues with their percentages
subset_eigenvals <- eigenvals[1:4, ]

# Plotting the first 4 eigenvectors as an inset using the computed percentages
eigen_inset <- subset_eigenvals %>%
  mutate(Component = row_number()) %>%
  ggplot(aes(x = Component, y = percentage)) +
  geom_bar(stat = "identity") +
  labs(title = "Eigenvalues",
       x = "Component",
       y = "Eigenvalue (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine the plots
combined <- grid.arrange(pca_plot, eigen_inset, ncol = 2)

# Save the plot
ggsave(filename = snakemake@output[['pca_plot']], plot = combined , width = 12, height = 6)
