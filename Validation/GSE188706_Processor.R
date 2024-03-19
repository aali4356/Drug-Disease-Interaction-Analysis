#Load the necessary libraries
library(GEOquery)
library(DESeq2)


# Opening and Unpacking Data ----------------------------------------------

# Set cwd as this folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("chdir.R")
source('nipals.R')
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

# Convert Ensembl IDs to Gene Symbols for GSE188706
symbols_GSE188706 <- mapIds(org.Hs.eg.db, 
                            keys=row.names(GSE188706.df.normCounts), 
                            column="SYMBOL", 
                            keytype="ENSEMBL", 
                            multiVals="first")
GSE188706.df.normCounts$symbols <- symbols_GSE188706


# Trim and Unique ---------------------------------------------------------

df_cleaned <- na.omit(GSE188706.df.normCounts, subset = "symbols")
# Get the vector of unique gene symbols
unique_gene_symbols <- unique(df_cleaned$symbols)

df_unique_symbols <- df_cleaned %>% distinct(symbols, .keep_all = TRUE)
rownames(df_unique_symbols) <- df_unique_symbols$symbols
df_unique_symbols <- df_unique_symbols[1:(length(df_unique_symbols)-2)]
df_unique_symbols <- df_unique_symbols[1:(length(df_unique_symbols)-3)]
new_column_names <- c(rep("0", 3), rep("1", 3))
names(df_unique_symbols) <- new_column_names

# Define the file path
file_path <- "GSE188706_Processed.txt"

# Use write.table to save the matrix to a file
write.table(df_unique_symbols, 
            file = file_path, sep = "\t", quote = FALSE, 
            row.names = TRUE)
