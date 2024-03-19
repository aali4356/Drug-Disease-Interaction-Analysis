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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Read the three text files and create data frames
file1 <- read_file("chdir_GSE188706_results.txt") # low dose G1
file2 <- read_file("chdir_GSE188706_High_results.txt") # high dose G1
file3 <- read_file("chdir_GSE130397_results.txt") # ER data

file1.ordered <- file1[order(file1$Value),]
file2.ordered <- file2[order(file2$Value),]
file3.ordered <- file3[order(file3$Value),]
n = 5000
file1 <- file1.ordered[c(1:n,(nrow(file1.ordered)-n):nrow(file1.ordered)),]
file2 <- file2.ordered[c(1:n,(nrow(file2.ordered)-n):nrow(file2.ordered)),]
file3 <- file3.ordered[c(1:n,(nrow(file3.ordered)-n):nrow(file3.ordered)),]
# Combine the data frames by the common 'Gene' column
combined_data <- merge(file1, file2, by = "Gene")
combined_data <- merge(combined_data, file3, by = "Gene")

# Replace NA values with 0
combined_data[is.na(combined_data)] <- 0

# Rename the columns
names(combined_data) <- c("Gene", "DElist1", "DElist2", "DElist3")

# Create the datatable
datatable <- combined_data[, c("DElist1", "DElist2", "DElist3")]

# Load the lsa library
library(lsa)

# Calculate the cosine similarities
cosine_sim_12 <- 1 - cosine(datatable[, 1], datatable[, 2])
cosine_sim_13 <- 1 - cosine(datatable[, 1], datatable[, 3])
cosine_sim_23 <- 1 - cosine(datatable[, 2], datatable[, 3])

# cosine_sim_12 <- 1 - cor(datatable[, 1], datatable[, 2], method = "spearman")
# cosine_sim_13 <- 1 - cor(datatable[, 1], datatable[, 3], method = "spearman")
# cosine_sim_23 <- 1 - cor(datatable[, 2], datatable[, 3], method = "spearman")


# Print the cosine similarities
cat("Cosine similarity between low G1 and high G1:", cosine_sim_12, "\n")
cat("Cosine similarity between low G1 and ER data:", cosine_sim_13, "\n")
cat("Cosine similarity between high G1 and ER data:", cosine_sim_23, "\n")
