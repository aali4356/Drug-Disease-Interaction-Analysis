# Extract the gene expression data and sample information
expression_data <- processed_glom_all$final
sample_info <- processed_glom_all$info$Group
platform_info <- processed_glom_all$info$Platform

# Filter out GPL570 data
gpl96_indices <- which(platform_info == "GPL96")
expression_data <- expression_data[, gpl96_indices]
sample_info <- sample_info[gpl96_indices]

# Get the unique disease groups
disease_groups <- unique(sample_info)

# Remove "Healthy Living Donor" and "Healthy Living Donor 570" from the disease groups
disease_groups <- disease_groups[!disease_groups %in% c("Healthy Living Donor", "Healthy Living Donor 570")]

# Combine "Healthy Living Donor" and "Healthy Living Donor 570" samples
healthy_indices <- which(sample_info %in% c("Healthy Living Donor", "Healthy Living Donor 570"))
healthy_data <- expression_data[, healthy_indices]

# Create a directory to store the text files
path <- "C:/Users/ahmad/OneDrive/BrandeisUniversity/__Research_120__/Data/Processed_Disease2"
dir.create(path, showWarnings = FALSE)

# Iterate over each disease group
for (disease in disease_groups) {
  # Get the indices of the current disease samples
  disease_indices <- which(sample_info == disease)
  
  # Combine healthy and disease data
  combined_data <- as.data.frame(cbind(healthy_data, expression_data[, disease_indices]))
  
  # Create labels for healthy (0) and disease (1) samples
  labels <- c(rep(0, ncol(healthy_data)), rep(1, length(disease_indices)))
  
  # Create a data frame with labels and expression data
  names(combined_data) <- labels
  combined_data <- na.omit(combined_data)
  
  # Generate the output file name
  output_file <- file.path(path, paste0(gsub("[^[:alnum:]]", "_", disease), "_vs_healthy.txt"))
  
  # Export the data frame as a text file
  write.table(combined_data, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE)
  
  cat("Exported:", output_file, "\n")
}