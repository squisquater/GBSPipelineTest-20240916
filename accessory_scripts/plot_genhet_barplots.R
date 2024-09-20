#!/usr/bin/env Rscript

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(viridis)  # Add the viridis package for color scaling

# Input and output file paths
input_file <- snakemake@input$merged
output_plot <- snakemake@output$plot
output_summary <- snakemake@output$summary  # Add an output for the summary .txt file

# Load the merged data
data <- read.table(input_file, header = TRUE, sep = "\t")

# List of genetic statistics to plot
metrics <- c("PHt", "Hs_obs", "Hs_exp", "IR", "HL")

# Function to create barplots of means and standard deviations
create_barplot <- function(data, metric) {
  data_summary <- data %>%
    group_by(Region.ID) %>%
    summarise(
      mean_value = mean(get(metric), na.rm = TRUE),
      sd_value = sd(get(metric), na.rm = TRUE),
      se_value = sd_value / sqrt(n())  # Calculate SE
    )
  
  ggplot(data_summary, aes(x = Region.ID, y = mean_value, fill = Region.ID)) +
    geom_bar(stat = "identity", color = "black") +
    geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = 0.2) +
    labs(x = "Population (Region.ID)", y = metric, title = paste("Comparison of", metric)) +
    theme_minimal() +
    scale_fill_hue(name = "Population") +  # Use scale_color_hue for discrete colors
    theme(legend.position = "none")
}

# Function to summarize the means and SE for each population and metric
summarize_metrics <- function(data, metrics) {
  summary_list <- lapply(metrics, function(metric) {
    data %>%
      group_by(Region.ID) %>%
      summarise(
        mean_value = mean(get(metric), na.rm = TRUE),
        sd_value = sd(get(metric), na.rm = TRUE),
        se_value = sd_value / sqrt(n())  # Calculate SE
      ) %>%
      mutate(metric = metric)
  })
  # Combine all summaries into a single data frame
  summary_df <- do.call(rbind, summary_list)
  return(summary_df)
}

# Create individual barplots for each metric
barplot_list <- lapply(metrics, function(metric) {
  create_barplot(data, metric)
})

# Combine all plots into a single figure using cowplot::plot_grid
combined_barplot <- plot_grid(plotlist = barplot_list, labels = "AUTO", ncol = 2)

# Save the combined plot as a single .png file
ggsave(filename = output_plot, plot = combined_barplot, width = 12, height = 10)

# Generate and save summary of means and SE for each metric
summary_data <- summarize_metrics(data, metrics)
write.table(summary_data, file = output_summary, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Combined barplot saved to:", output_plot, "\n")
cat("Summary of means and SE saved to:", output_summary, "\n")
