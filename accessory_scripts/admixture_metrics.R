# admixture_metrics.R

input_dir <- snakemake@params$input_dir
input_prefix <- snakemake@params$input_prefix
output_file <- snakemake@output[[1]]
K_start <- as.integer(snakemake@params$K_start)
K_end <- as.integer(snakemake@params$K_end)

# Print out the input directory and prefix for debugging
cat("Input directory:", input_dir, "\n")
cat("Input prefix:", input_prefix, "\n")

results <- data.frame(K = integer(), Loglikelihood = numeric(), CV_error = numeric(), stringsAsFactors = FALSE)

for (K in K_start:K_end) {
  # Construct the file path
  file_path <- sprintf("%s/%s_renamedChr_K%d.Q", input_dir, input_prefix, K)
  
  # Print the file path to check if it's correct
  cat("Attempting to read file:", file_path, "\n")
  
  # Try reading the file
  if (file.exists(file_path)) {
    file_content <- readLines(file_path)
    
    # Check if the file content was read successfully
    if (length(file_content) == 0) {
      cat("Warning: The file", file_path, "is empty.\n")
    }
    
    summary_index <- which(grepl("Summary:", file_content))
    
    if (length(summary_index) > 0) {
      # Extract the log-likelihood
      loglikelihood_line <- file_content[summary_index + 2] # Assuming Loglikelihood is two lines below "Summary:"
      loglikelihood <- as.numeric(gsub("Loglikelihood: ([0-9.-]+)", "\\1", loglikelihood_line))
      
      # Find and extract the CV error
      cv_error <- NA # Default in case CV error is not found
      for (i in (summary_index + 3):length(file_content)) {
        if (grepl("CV error \\(K=", file_content[i])) {
          cv_error_line <- file_content[i]
          cv_error <- as.numeric(gsub(".*CV error \\(K=\\d+\\): ([0-9.]+)", "\\1", cv_error_line))
          break # Stop the loop once the CV error line is found
        }
      }
    } else {
      loglikelihood <- NA
      cv_error <- NA
      cat("Warning: No 'Summary:' found in file", file_path, "\n")
    }
  } else {
    cat("Error: File", file_path, "does not exist.\n")
    loglikelihood <- NA
    cv_error <- NA
  }
  
  results <- rbind(results, data.frame(K = K, Loglikelihood = loglikelihood, CV_error = cv_error))
}

# Write the results to the output file
write.table(results, file = output_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
