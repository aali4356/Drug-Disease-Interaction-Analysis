#Load the necessary libraries
library(GEOquery)
library(DESeq2)


# Opening and Unpacking Data ----------------------------------------------

# Set cwd as this folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("chdir.R")
library(readr)
library(AnnotationDbi)
library(hgu133plus2.db)
GSE188706.file.data <-  "GSE188706_Raw_gene_counts_matrix.csv"
data.raw <- read_csv(GSE188706.file.data)
GSE188706.df <- as.data.frame(data.raw)


# Normalization -----------------------------------------------------------


colData <- data.frame(
  sampleName = colnames(GSE188706.df)[-1], # Assuming first column is gene IDs
  condition = c(rep("control", 3), rep("low_dose", 3), rep("high_dose", 3))
)
rownames(colData) <- colData$sampleName

# Make sure gene names (Ensembl IDs) are the row names of the count data
rownames(GSE188706.df) <- GSE188706.df[,1]
# Remove the gene names column before creating DESeqDataSet
countData <- GSE188706.df[,-1]

# Create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)

# Run the DESeq normalization and differential expression
dds <- DESeq(dds,)

# Getting normalized counts
normCounts <- counts(dds, normalized = TRUE)

# Conducting a differential expression analysis between two conditions, for example
GSE188706.df.results <- results(dds, contrast=c("condition", "low_dose", "control"))
# Add gene names as a column to the normCounts matrix for reference or further processing
new_column_names <- c(rep("0", ), rep("1", 4))
GSE188706.df.normCounts <- as.data.frame(normCounts)
GSE188706.df.normCounts$Ensembl_ID <- rownames(GSE188706.df.normCounts)

highdose.g1 <- GSE188706.df.normCounts[,c(1,2,3,7,8,9)]
new_column_names <- c(rep("0", 3), rep("1",3))
names(highdose.g1) <- new_column_names
# Second set of data ------------------------------------------------------

library(dplyr)
library(readr)
library(purrr)
# Define the paths to your files within the "folder"
files_experimental <- list.files(path = "GSE130397_RAW", pattern = "GSM3737480|GSM3737479|GSM3737478|GSM3737477|GSM3737475|GSM3737473", full.names = TRUE)
files_control <- list.files(path = "GSE130397_RAW", pattern = "GSM3737472|GSM3737471|GSM3737470|GSM3737469|GSM3737468|GSM3737467", full.names = TRUE)

# Function to process each file
process_file <- function(file_path) {
  data <- read_tsv(file_path, col_names = FALSE, col_types = cols(.default = "c")) %>% 
    select(X1, X2) %>% 
    rename(GeneID = X1, Count = X2)
  
  # Extract sample name from file name
  sample_name <- gsub(".*_(S\\d+_\\d+).*", "\\1", basename(file_path))
  colnames(data)[2] <- sample_name
  
  return(data)
}

# Process experimental files
data_experimental <- lapply(files_experimental, process_file) %>% 
  purrr::reduce(full_join, by = "GeneID")

# Process control files
data_control <- lapply(files_control, process_file) %>% 
  purrr::reduce(full_join, by = "GeneID")

GSE130397.data <- cbind(data_control, data_experimental[,-1])
GSE130397.data <- GSE130397.data[-1,]
rownames(GSE130397.data) <- GSE130397.data[,1]
GSE130397.data <- as.data.frame(GSE130397.data[,-1])
colnames(GSE130397.data) <- c("S1", "S2", "S3", "S4", "S5", "S6",
                              "ER1", "ER2", "ER3", "ER4", "ER5", "ER6")
# Second Normalization ----------------------------------------------------


colData.2 <- data.frame(
  sampleName = colnames(GSE130397.data), 
  condition = c(rep("control", 6), rep("experimental", 6))
)
rownames(colData.2) <- colData.2$sampleName

# Remove the gene names column before creating DESeqDataSet

mat_fun <- function(m){
  m2 <- apply(m,  2,  function(x) as.numeric(paste(x)))
  colnames(m2) <- colnames(m)
  rownames(m2) <- rownames(m)
  return(m2)
}

countData <- mat_fun(GSE130397.data)
# Run the DESeq normalization and differential expression

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData.2,
                              design = ~ condition)
dds <- DESeq(dds)

# Getting normalized counts
normCounts <- counts(dds, normalized = TRUE)

# Conducting a differential expression analysis between two conditions, for example
results <- results(dds, contrast=c("condition", "experimental", "control"))
# Add gene names as a column to the normCounts matrix for reference or further processing
GSE130397.normCounts <- as.data.frame(normCounts)
rownames(GSE130397.normCounts) <- rownames(GSE130397.data)
GSE130397.normCounts <- as.data.frame(GSE130397.normCounts)
new_column_names <- c(rep("0", 6), rep("1", 6))
names(GSE130397.normCounts) <- new_column_names
# Define the file path
file_path <- "GSE130397_Processed.txt"

# Use write.table to save the matrix to a file
write.table(GSE130397.normCounts, 
            file = file_path, sep = "\t", quote = FALSE, 
            row.names = TRUE)

# GSE188706.df.normCounts
GSE188706.df.normCounts <- as.data.frame(GSE188706.df.normCounts[,c(1,2,3,4,5,6)])
new_column_names <- c(rep("0", 3), rep("1", 3))
names(GSE188706.df.normCounts) <- new_column_names
file_path <- "GSE188706_Processed_ensembl.txt"
write.table(GSE188706.df.normCounts, 
            file = file_path, sep = "\t", quote = FALSE, 
            row.names = TRUE)


# chdir(GSE188706.df.normCounts)
# Gene Conversions --------------------------------------------------------

# Convert Ensembl IDs to Gene Symbols for GSE188706
symbols_GSE188706 <- mapIds(org.Hs.eg.db, 
                            keys=row.names(GSE188706.df.normCounts), 
                            column="SYMBOL", 
                            keytype="ENSEMBL", 
                            multiVals="first")

# High gene G1 Conversion
symbols_GSE188706_high <- mapIds(org.Hs.eg.db, 
                            keys=row.names(highdose.g1), 
                            column="SYMBOL", 
                            keytype="ENSEMBL", 
                            multiVals="first")
# Replace row names with Gene Symbols
# row.names(highdose.g1) <- symbols_GSE188706_high
# Do the same for GSE130397
symbols_GSE130397 <- mapIds(org.Hs.eg.db, 
                            keys=row.names(GSE130397.normCounts), 
                            column="SYMBOL", 
                            keytype="ENSEMBL", 
                            multiVals="first")

# Replace row names with Gene Symbols
# row.names(GSE130397.normCounts) <- symbols_GSE130397
# Add the gene symbols directly as a new column
GSE188706.df.normCounts$symbols <- symbols_GSE188706
GSE130397.normCounts$symbols <- symbols_GSE130397
highdose.g1$symbols <- symbols_GSE188706_high


# Function to append a suffix to duplicate values to make them unique
make_unique <- function(x) {
  ux <- unique(x)
  for(i in ux) {
    wh <- which(x == i)
    if(length(wh) > 1) x[wh] <- paste(i, seq_along(wh), sep="_")
  }
  x
}

# Apply this function to your symbols vector before setting row names
unique_symbols_GSE188706 <- make_unique(as.character(symbols_GSE188706))
unique_symbols_GSE188706.high <- make_unique(as.character(symbols_GSE188706_high))

unique_symbols_GSE130397 <- make_unique(as.character(symbols_GSE130397))
names(GSE130397.normCounts) <- c("S1", "S2", "S3", "S4", "S5", "S6",
                                 "ER1", "ER2", "ER3", "ER4", "ER5", "ER6", "symbols")
# Create a temporary DataFrame with row names (Ensembl IDs) and symbols
temp_df <- data.frame(Ensembl_ID = row.names(GSE130397.normCounts), Symbol = symbols_GSE130397, stringsAsFactors = FALSE)

# Bind the symbol information to the original count data
GSE130397.normCountss <- as.data.frame(cbind(temp_df, GSE130397.normCounts))

# Aggregate counts by gene symbol, taking the mean (or another appropriate summary statistic) of duplicate symbols
GSE130397.normCounts.aggregated <- GSE130397.normCountss %>%
  group_by(Symbol) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
library(dplyr)

GSE130397.normCounts.cleaned <- GSE130397.normCountss %>%
  filter(!is.na(Symbol))
# GSE130397.normCounts.aggregated <- GSE130397.normCounts.cleaned %>%
#   group_by(Symbol) %>%
#   summarise(across(everything(), ~mean(.x, na.rm = TRUE)))
backup <- GSE130397.normCounts.cleaned
GSE130397.normCounts.cleaned <- GSE130397.normCounts.cleaned %>%
  distinct(Symbol, .keep_all = TRUE)


# If you then want to set these symbols as row names
row.names(GSE130397.normCounts.cleaned) <- GSE130397.normCounts.cleaned$Symbol

GSE130397.normCounts.cleanedd <- GSE130397.normCounts.cleaned[,c(3:14)]
names(GSE130397.normCounts.cleanedd) <- new_column_names
file_path <- "GSE130397_Processed.txt"

# Use write.table to save the matrix to a file
write.table(GSE130397.normCounts.cleanedd, 
            file = file_path, sep = "\t", quote = FALSE, 
            row.names = TRUE)

file_path <- "GSE188706_High_Processed.txt"
highdose.g1 <- na.omit(highdose.g1, subset = symbols)
# highdose.g1 <- highdose.g1 %>% filter(!is.na(symbols))
highdose.g1 <- highdose.g1[!duplicated(highdose.g1$symbols), ]
rownames(highdose.g1) <- highdose.g1$symbols
highdose.g1$symbols <- NULL

# Use write.table to save the matrix to a file
write.table(highdose.g1, 
            file = file_path, sep = "\t", quote = FALSE, 
            row.names = TRUE)
