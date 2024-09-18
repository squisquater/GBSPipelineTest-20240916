#!/usr/bin/env Rscript

# Load necessary libraries
library(ggplot2)
library(cowplot)
library(dplyr)

# Input and output file paths
merged_file <- snakemake@input$merged
plot_output <- snakemake@output$plot  # Single PNG file to save the combined plot
results_output <- snakemake@output$results  # Text file to save correlation results

# Load the merged data
merged_data <- read.table(merged_file, header = TRUE, sep = "\t")

# Check if merged data loaded successfully
cat("Merged data loaded successfully\n")

# List of metrics to plot
metrics <- c("PHt", "Hs_obs", "Hs_exp", "IR", "HL")

# Create a function to run correlation analysis and return statistics
calculate_stats <- function(df, x_var, y_var) {
  # Perform linear regression
  model <- lm(get(y_var) ~ get(x_var), data = df)
  summary_model <- summary(model)
  
  # Extract the slope and p-value
  slope <- summary_model$coefficients[2, 1]
  p_value <- summary_model$coefficients[2, 4]
  r_value <- cor(df[[x_var]], df[[y_var]], use = "complete.obs")
  
  return(data.frame(metric = y_var, Region.ID = unique(df$Region.ID), slope = slope, r_value = r_value, p_value = p_value))
}

# Create individual plots for each metric with color coding for Region.ID (population)
plot_list <- lapply(metrics, function(metric) {
  ggplot(merged_data, aes(x = GenomeWideAvgDepth, y = get(metric), color = Region.ID)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "Genome-Wide Average Depth", y = metric, title = paste("Depth vs", metric)) +
    theme_minimal() +
    scale_color_discrete(name = "Population")
})

# Combine all plots into a single figure using cowplot::plot_grid
combined_plot <- plot_grid(plotlist = plot_list, labels = "AUTO", ncol = 2)

# Save the combined plot as a single .png file
ggsave(filename = plot_output, plot = combined_plot, width = 10, height = 10)

# Run correlation analysis for each metric within each population
stats_results <- merged_data %>%
  group_by(Region.ID) %>%
  do({
    results <- lapply(metrics, function(metric) {
      calculate_stats(., x_var = "GenomeWideAvgDepth", y_var = metric)
    })
    bind_rows(results)
  })

# Print the correlation results to check
print(stats_results)

# Check if stats_results is not empty
if (nrow(stats_results) > 0) {
  # Save the results to a separate text file
  write.table(stats_results, file = results_output, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Correlation statistics saved to:", results_output, "\n")
} else {
  cat("No correlation statistics calculated. Please check the data.\n")
}
