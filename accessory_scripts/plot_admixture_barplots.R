# Load necessary libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(tidyr)

# Read arguments from the snakemake object
admixture_dir <- snakemake@params$admixture_dir
summary_file <- snakemake@input$summary_file
output_file <- snakemake@output$output_file

# Load the summary data
data <- read.table(summary_file, header = TRUE, sep = "\t")

# Gather data into long format for ggplot
data_long <- data %>%
  pivot_longer(cols = starts_with("K"), names_to = "Cluster", values_to = "Proportion") %>%
  mutate(K_value = as.numeric(gsub("K(\\d+)_.*", "\\1", Cluster)))

# Initialize a list to store plots
plots <- list()

# Loop over unique K values and create a plot for each
for (K in sort(unique(data_long$K_value))) {
  # Filter data for the current K value
  data_K <- data_long %>% filter(K_value == K)
  
  # Create a stacked bar plot for this K value
  p <- ggplot(data_K, aes(x = Individual_ID, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste("K =", K), x = "Individual", y = "Proportion") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
  
  # Add the plot to the list
  plots[[K]] <- p
}

# Combine all the plots into a single plot using cowplot
final_plot <- plot_grid(plotlist = plots, ncol = 1)

# Save the combined plot to a file
ggsave(output_file, final_plot, width = 10, height = length(plots) * 4)
