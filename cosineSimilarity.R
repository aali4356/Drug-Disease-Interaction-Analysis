# Load the required libraries
library(dplyr)
library(tidyr)
library(lsa)
library(pheatmap)

# Function to read the text files and create a data frame
read_file <- function(file_path) {
  lines <- readLines(file_path)
  gene_names <- character()
  values <- numeric()
  for (line in lines) {
    parts <- strsplit(line, ", ")[[1]]
    gene_names <- c(gene_names, parts[1])
    values <- c(values, as.numeric(parts[2]))
  }
  data.frame(Gene = gene_names, Value = values, stringsAsFactors = FALSE)
}

# Set the working directory to the location of the script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Get the list of drug text files in the "Data/Data_Output" folder
drug_file_list <- list.files("Data/Data_Output", pattern = "*.txt", full.names = TRUE)

# Get the list of disease text files in the "Data/Disease_Data_Output" folder
disease_file_list <- list.files("Data/Disease_Data_Output", pattern = "*.txt", full.names = TRUE)

# Initialize empty lists to store the data frames
drug_data_list <- list()
disease_data_list <- list()

# Read each drug text file and create a data frame
for (file_path in drug_file_list) {
  file_name <- gsub("_results\\.txt$", "", basename(file_path))
  data_frame <- read_file(file_path)
  data_frame$File <- file_name
  drug_data_list[[file_name]] <- data_frame
}

# Read each disease text file and create a data frame
for (file_path in disease_file_list) {
  file_name <- gsub("chdir_|_vs_healthy_results\\.txt$", "", basename(file_path))
  data_frame <- read_file(file_path)
  data_frame$File <- file_name
  disease_data_list[[file_name]] <- data_frame
}

# Combine all drug data frames into a single data frame
combined_drug_data <- bind_rows(drug_data_list, .id = "File")

# Combine all disease data frames into a single data frame
combined_disease_data <- bind_rows(disease_data_list, .id = "File")

# Convert gene names to uppercase
combined_drug_data$Gene <- toupper(combined_drug_data$Gene)
combined_disease_data$Gene <- toupper(combined_disease_data$Gene)

# Find the common genes across all files (drug and disease)
common_genes <- intersect(combined_drug_data$Gene, combined_disease_data$Gene)

# Filter the combined data to keep only the common genes
filtered_drug_data <- combined_drug_data %>%
  filter(Gene %in% common_genes) %>%
  pivot_wider(names_from = File, values_from = Value)

filtered_disease_data <- combined_disease_data %>%
  filter(Gene %in% common_genes) %>%
  pivot_wider(names_from = File, values_from = Value)

# Create the data matrix for cosine similarity calculation
data_matrix <- as.matrix(filtered_drug_data[, -1])
rownames(data_matrix) <- filtered_drug_data$Gene

# Add the disease data to the data matrix
data_matrix <- cbind(data_matrix, as.matrix(filtered_disease_data[, -1]))
data_matrix <- na.omit(data_matrix)
# Calculate the cosine similarity matrix
cosine_sim_matrix <- matrix(0, nrow = ncol(data_matrix), ncol = ncol(data_matrix))
for (i in 1:(ncol(data_matrix) - 1)) {
  for (j in (i + 1):ncol(data_matrix)) {
    cosine_sim_matrix[i, j] <- 1 - cosine(data_matrix[, i], data_matrix[, j])
    cosine_sim_matrix[j, i] <- cosine_sim_matrix[i, j]  # Symmetric matrix
  }
}

# Set the column and row names of the cosine similarity matrix
colnames(cosine_sim_matrix) <- colnames(data_matrix)
rownames(cosine_sim_matrix) <- colnames(data_matrix)

subset_matrix <- cosine_sim_matrix[c(1:24),c(25:32)]
subset_matrix_cleaned <- subset_matrix[-c(12,18,20,24),]
drug_matrix_cleaned <- cosine_sim_matrix[c(1:11,13,15:24), c(1:11,13,15:24)]
pheatmap(subset_matrix_cleaned, main = "Cosine Similarity Heatmap", fontsize = 8, angle_col = 45)


# Create a heatmap of the cosine similarity matrix
pheatmap(cosine_sim_matrix, main = "Cosine Similarity Heatmap", fontsize = 8, angle_col = 45)



# Further analysis --------------------------------------------------------
# Calculate the cosine similarity matrix
library(lsa)
cosine_sim_matrix <- matrix(0, nrow = ncol(data_matrix), ncol = ncol(data_matrix))
for (i in 1:(ncol(data_matrix) - 1)) {
  for (j in (i + 1):ncol(data_matrix)) {
    cosine_sim_matrix[i, j] <- 1 - cosine(data_matrix[, i], data_matrix[, j])
    cosine_sim_matrix[j, i] <- cosine_sim_matrix[i, j] # Symmetric matrix
  }
}

# Set the column and row names of the cosine similarity matrix
colnames(cosine_sim_matrix) <- colnames(data_matrix)
rownames(cosine_sim_matrix) <- colnames(data_matrix)

subset_matrix <- cosine_sim_matrix[c(1:24),c(25:32)]
subset_matrix_cleaned <- subset_matrix[-c(12,18,20,24),]

# Perform permutation test on the subset_matrix_cleaned
num_permutations <- 1000
p_values <- matrix(0, nrow = nrow(subset_matrix_cleaned), ncol = ncol(subset_matrix_cleaned))

for (i in 1:nrow(subset_matrix_cleaned)) {
  for (j in 1:ncol(subset_matrix_cleaned)) {
    observed_cosine_sim <- subset_matrix_cleaned[i, j]
    
    # Perform permutations
    permuted_cosine_sims <- rep(0, num_permutations)
    for (k in 1:160) {
      permuted_values <- sample(subset_matrix_cleaned)
      permuted_cosine_sim <- permuted_values[k]
      permuted_cosine_sims[k] <- permuted_cosine_sim
    }
    
    # Calculate the p-value
    p_values[i, j] <- sum(permuted_cosine_sims >= observed_cosine_sim) / num_permutations
  }
}

# Create a matrix of asterisks based on the p-values
significance_matrix <- matrix("", nrow = nrow(p_values), ncol = ncol(p_values))
significance_matrix[p_values < 0.05] <- "*"
significance_matrix[p_values < 0.01] <- "**"
significance_matrix[p_values < 0.001] <- "***"

# Create the heatmap with significance asterisks
dev.off()
heatmap.2(subset_matrix_cleaned,
          main = "Cosine Similarity Heatmap",
          cexRow = 0.8,
          cexCol = 0.8,
          key = TRUE,
          symkey = FALSE,
          density.info = "none",
          trace = "none",
          cellnote = significance_matrix,
          notecol = "black",
          notecex = 1.4,
          col = colorRampPalette(c("blue", "white", "red"))(100),
          margins=c(11,14)
          )

# Perform permutation test on the subset_matrix_cleaned
num_permutations <- 1000
p_values <- matrix(0, nrow = nrow(drug_matrix_cleaned), ncol = ncol(drug_matrix_cleaned))

for (i in 1:nrow(drug_matrix_cleaned)) {
  for (j in 1:ncol(drug_matrix_cleaned)) {
    observed_cosine_sim <- drug_matrix_cleaned[i, j]
    
    # Perform permutations
    permuted_cosine_sims <- rep(0, num_permutations)
    for (k in 1:484) {
      permuted_values <- sample(drug_matrix_cleaned)
      permuted_cosine_sim <- permuted_values[k]
      permuted_cosine_sims[k] <- permuted_cosine_sim
    }
    
    # Calculate the p-value
    p_values[i, j] <- sum(permuted_cosine_sims >= observed_cosine_sim) / num_permutations
  }
}

# Create a matrix of asterisks based on the p-values
significance_matrix <- matrix("", nrow = nrow(p_values), ncol = ncol(p_values))
significance_matrix[p_values < 0.05] <- "*"
significance_matrix[p_values < 0.01] <- "**"
significance_matrix[p_values < 0.001] <- "***"

heatmap.2(drug_matrix_cleaned,
          main = "Cosine Similarity Heatmap",
          # cexRow = 0.8,
          # cexCol = 0.8,
          key = TRUE,
          symkey = FALSE,
          # srtCol=45, adjCol = c(1,1),
          # srtRow=45, adjRow=c(0, 1),
          density.info = "none",
          trace = "none",
          cellnote = significance_matrix,
          notecol = "black",
          notecex = 1.4,
          col = colorRampPalette(c("blue", "white", "red"))(100),
          margins=c(13,14)
)
