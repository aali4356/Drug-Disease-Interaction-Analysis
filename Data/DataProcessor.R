#Load the necessary libraries
library(GEOquery)
library(DESeq2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# source("chdir.R")
# source('nipals.R')
library(readr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(dplyr)


# Atorvastatin ------------------------------------------------------------

# Opening and Unpacking Data


atorvastatin.file.data <-  "Atorvastatin_GSE196701/GSE196701_gene_expression.csv"
atorvastatin.data.raw <- read_csv(atorvastatin.file.data)
atorvastatin.data.raw <- as.data.frame(atorvastatin.data.raw)


# Cleaning Data

atorvastatin.data <- atorvastatin.data.raw[,c(3:9)]
atorvastatin.data.cleaned <- na.omit(atorvastatin.data, subset = "gene_symbol")

atorvastatin.data.cleaned <- atorvastatin.data.cleaned %>%
  distinct(gene_symbol, .keep_all = TRUE)
row.names(atorvastatin.data.cleaned) <- atorvastatin.data.cleaned$gene_symbol

ato.df <- cbind(atorvastatin.data.cleaned[,c(5:7)], 
                atorvastatin.data.cleaned[,c(2:4)])
Atorvastatin_GSE196701_Processed <- ato.df
new_column_names <- c(rep("0", 3), rep("1", 3))
names(Atorvastatin_GSE196701_Processed) <- new_column_names

# Write on File 

# Define the file path
file_path <- "Processed/Atorvastatin_GSE196701_Processed.txt"

# Use write.table to save the matrix to a file
write.table(Atorvastatin_GSE196701_Processed, 
            file = file_path, sep = "\t", quote = FALSE, 
            row.names = TRUE)

# Bicarbonate  ------------------------------------------------------------

# Opening and Unpacking Data

library(org.Mm.eg.db)

bicarbonate.file.data <-  "Bicarbonate_GSE200638/GSE200638_rpkm_IRR.csv.gz"
bicarbonate.data.raw <- read_csv(bicarbonate.file.data)
bicarbonate.data.raw <- as.data.frame(bicarbonate.data.raw)


# Cleaning Data

bicarbonate.data <- bicarbonate.data.raw[,c(1:9)]
# bicarbonate.data.cleaned <- na.omit(bicarbonate.data, subset = "gene_symbol")
names(bicarbonate.data) <- c("gene_id", 
                             "wt1", "wt2", "wt3", "wt4",
                             "bc1", "bc2", "bc3", "bc4") 
bicarbonate.data.cleaned <- bicarbonate.data %>%
  distinct(gene_id, .keep_all = TRUE)

gene_symbols <- AnnotationDbi::select(org.Mm.eg.db, keys = bicarbonate.data.cleaned$gene_id, 
                       columns = c("SYMBOL"), 
                       keytype = "ENSEMBL")

gene_symbols <- gene_symbols %>%
  distinct(SYMBOL, .keep_all = TRUE)

names(bicarbonate.data.cleaned)[names(bicarbonate.data.cleaned) == "gene_id"] <- "ENSEMBL"

bicarbonate.data.cleaned_with_symbols <- 
  merge(bicarbonate.data.cleaned, gene_symbols, by = "ENSEMBL", all.x = TRUE, no.dups = TRUE)

bicarbonate.data.cleaned <- bicarbonate.data.cleaned_with_symbols
bicarbonate.data.cleaned_unique <- bicarbonate.data.cleaned[!duplicated(bicarbonate.data.cleaned$SYMBOL), ]
bicarbonate.data.cleaned_unique <- bicarbonate.data.cleaned_unique[!is.na(bicarbonate.data.cleaned_unique$SYMBOL), ]
row.names(bicarbonate.data.cleaned_unique) <- bicarbonate.data.cleaned_unique$SYMBOL

bicarbonate.data.cleaned_unique$SYMBOL <- NULL
bicarbonate.data.cleaned_unique$ENSEMBL <- NULL
bicarb.df <- bicarbonate.data.cleaned_unique



Bicarbonate_GSE200638_Processed <- bicarb.df
new_column_names <- c(rep("0", 4), rep("1", 4))
names(Bicarbonate_GSE200638_Processed) <- new_column_names

# Write on File 

# Define the file path
file_path <- "Processed/Bicarbonate_GSE200638_Processed.txt"

# Use write.table to save the matrix to a file
write.table(Bicarbonate_GSE200638_Processed, 
            file = file_path, sep = "\t", quote = FALSE, 
            row.names = TRUE)



# Fineronone  ------------------------------------------------------------

# Opening and Unpacking Data

fineronone.file.data <-  "Fineronone_GSE183841/GSE183841_EXPORT_BulkRNAseq_TPM_Counts.csv"
fineronone.data.raw <- read_csv(fineronone.file.data)
fineronone.data.raw <- as.data.frame(fineronone.data.raw)


# Cleaning Data

fineronone.data <- fineronone.data.raw
fineronone.data <- na.omit(fineronone.data, subset = "gene")
# names(fineronone.data) <- c("gene_id", 
#                              "wt1", "wt2", "wt3", "wt4",
#                              "bc1", "bc2", "bc3", "bc4") 
fineronone.data.cleaned <- fineronone.data %>%
  distinct(gene, .keep_all = TRUE)
library(org.Rn.eg.db)
gene_symbols <- AnnotationDbi::select(org.Rn.eg.db, keys = fineronone.data.cleaned$gene, 
                                      columns = c("SYMBOL"), 
                                      keytype = "ENSEMBL")

gene_symbols <- gene_symbols %>%
  distinct(SYMBOL, .keep_all = TRUE)

names(fineronone.data.cleaned)[names(fineronone.data.cleaned) == "gene"] <- "ENSEMBL"

fineronone.data.cleaned_with_symbols <- 
  merge(fineronone.data.cleaned, gene_symbols, by = "ENSEMBL", all.x = TRUE, no.dups = TRUE)

fineronone.data.cleaned <- fineronone.data.cleaned_with_symbols
fineronone.data.cleaned_unique <- fineronone.data.cleaned[!duplicated(fineronone.data.cleaned$SYMBOL), ]
fineronone.data.cleaned_unique <- fineronone.data.cleaned_unique[!is.na(fineronone.data.cleaned_unique$SYMBOL), ]
row.names(fineronone.data.cleaned_unique) <- fineronone.data.cleaned_unique$SYMBOL

fineronone.data.cleaned_unique$SYMBOL <- NULL
fineronone.data.cleaned_unique$ENSEMBL <- NULL
fineronone.df <- fineronone.data.cleaned_unique[,c(5:12)]



Fineronone_GSE183841_Processed <- fineronone.df
new_column_names <- c(rep("0", 4), rep("1", 4))
names(Fineronone_GSE183841_Processed) <- new_column_names

# Write on File 

# Define the file path
file_path <- "Processed/Fineronone_GSE183841_Processed.txt"

# Use write.table to save the matrix to a file
write.table(Fineronone_GSE183841_Processed, 
            file = file_path, sep = "\t", quote = FALSE, 
            row.names = TRUE)






# Losartan ----------------------------------------------------------------


# Opening and Unpacking Data


losartan.file.data <-  "Losartan_ARB_GSE159059/GSE159059_sorted_gene_fpkm_replaced.csv"
losartan.data.raw <- read_csv(losartan.file.data)
losartan.data.raw <- as.data.frame(losartan.data.raw)


# Cleaning Data

losartan.data <- losartan.data.raw[,c(1:21)]
losartan.data.cleaned <- na.omit(losartan.data, subset = "tracking_id")

losartan.data.cleaned <- losartan.data.cleaned %>%
  distinct(tracking_id, .keep_all = TRUE)
row.names(losartan.data.cleaned) <- losartan.data.cleaned$tracking_id
losartan.data.cleaned$tracking_id <- NULL
# los.df <- cbind(losartan.data.cleaned[,c(5:7)], 
#                 losartan.data.cleaned[,c(2:4)])
losartan_GSE159059_Processed <- losartan.data.cleaned
new_column_names <- c(rep("0", 10), rep("1", 10))
names(losartan_GSE159059_Processed) <- new_column_names

# Write on File 

# Define the file path
file_path <- "Processed/losartan_GSE159059_Processed.txt"

# Use write.table to save the matrix to a file
write.table(losartan_GSE159059_Processed, 
            file = file_path, sep = "\t", quote = FALSE, 
            row.names = TRUE)


# Mycophenolate Mofetil ---------------------------------------------------

# Opening and Unpacking Data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("MMF_GSE153021")
# mmf.file.folder <-  "MMF_GSE1563/GSE1563_RAW"
# mmf.data.raw <- read_csv(mmf.file.folder)
files <- list.files(pattern = "\\.txt$")
read_file <- function(file) {
  df <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  colnames(df)[2] <- sub("\\.txt$", "", basename(file)) # Optional: Rename the second column to the file name (or any other unique identifier)
  return(df)
}
# Read and merge all files
all_data <- lapply(files, read_file) %>% 
  Reduce(function(x, y) full_join(x, y, by = "Symbol"), .)

vehicle <- all_data[,c(4,6,7)]
mmf <- all_data[,c(2,3,5)]
genes <- all_data$Symbol
rebuilt.mmf.data <- data.frame(genes,
                              vehicle, 
                              mmf)
rebuilt.mmf.data <- rebuilt.mmf.data[!duplicated(rebuilt.mmf.data$genes), ]

rownames(rebuilt.mmf.data) <- rebuilt.mmf.data$genes
rebuilt.mmf.data$genes <- NULL
new_column_names <- c(rep("0", 3), rep("1", 3))
names(rebuilt.mmf.data) <- new_column_names

# Define the file path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
file_path <- "Processed/MMF_GSE153021_Processed.txt"

# Use write.table to save the matrix to a file
write.table(rebuilt.mmf.data, 
            file = file_path, sep = "\t", quote = FALSE, 
            row.names = TRUE)



# Prednisone --------------------------------------------------------------

# Opening and Unpacking Data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("Prednisone_GSE153021")
# mmf.file.folder <-  "MMF_GSE1563/GSE1563_RAW"
# mmf.data.raw <- read_csv(mmf.file.folder)
files <- list.files(pattern = "\\.txt$")
read_file <- function(file) {
  df <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  colnames(df)[2] <- sub("\\.txt$", "", basename(file)) # Optional: Rename the second column to the file name (or any other unique identifier)
  return(df)
}
# Read and merge all files
all_data <- lapply(files, read_file) %>% 
  Reduce(function(x, y) full_join(x, y, by = "Symbol"), .)

vehicle <- all_data[,c(2,6,7)]
pred <- all_data[,c(3,4,5)]
genes <- all_data$Symbol
rebuilt.pred.data <- data.frame(genes,
                               vehicle, 
                               pred)
rebuilt.pred.data <- rebuilt.pred.data[!duplicated(rebuilt.pred.data$genes), ]

rownames(rebuilt.pred.data) <- rebuilt.pred.data$genes
rebuilt.pred.data$genes <- NULL
new_column_names <- c(rep("0", 3), rep("1", 3))
names(rebuilt.pred.data) <- new_column_names

# Define the file path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
file_path <- "Processed/Prednisone_GSE153021_Processed.txt"

# Use write.table to save the matrix to a file
write.table(rebuilt.pred.data, 
            file = file_path, sep = "\t", quote = FALSE, 
            row.names = TRUE)


# Ramipril ----------------------------------------------------------------

# Ramipril_RASi_GSE226353
# Opening and Unpacking Data


ramipril.file.data <-  "Ramipril_RASi_GSE226353/GSE226353_ex_rpmk.csv.gz"
ramipril.data.raw <- read_csv(ramipril.file.data)
ramipril.data.raw <- as.data.frame(ramipril.data.raw)


# Cleaning Data

ramipril.data <- ramipril.data.raw[,c(2,7:12)]
ramipril.data.cleaned <- na.omit(ramipril.data, subset = "Geneid")

ramipril.data.cleaned <- ramipril.data.cleaned %>%
  distinct(Geneid, .keep_all = TRUE)
row.names(ramipril.data.cleaned) <- ramipril.data.cleaned$Geneid
ramipril.data.cleaned$Geneid <- NULL
# rami.df <- cbind(ramipril.data.cleaned[,c(5:7)], 
#                 ramipril.data.cleaned[,c(2:4)])
ramipril_GSE226353_Processed <- ramipril.data.cleaned
new_column_names <- c(rep("0", 3), rep("1", 3))
names(ramipril_GSE226353_Processed) <- new_column_names

# Write on File 

# Define the file path
file_path <- "Processed/ramipril_GSE226353_Processed.txt"

# Use write.table to save the matrix to a file
write.table(ramipril_GSE226353_Processed, 
            file = file_path, sep = "\t", quote = FALSE, 
            row.names = TRUE)

# Cleaning Ramipril -------------------------------------------------------
# Ramipril ----------------------------------------------------------------
# Ramipril_RASi_GSE226353

# Opening and Unpacking Data
ramipril.file.data <- "Ramipril_RASi_GSE226353/GSE226353_ex_rpmk.csv.gz"
ramipril.data.raw <- read_csv(ramipril.file.data)
ramipril.data.raw <- as.data.frame(ramipril.data.raw)

# Cleaning Data
ramipril.data <- ramipril.data.raw[,c(2,7:12)]
ramipril.data.cleaned <- na.omit(ramipril.data, subset = "Geneid")
ramipril.data.cleaned <- ramipril.data.cleaned %>%
  distinct(Geneid, .keep_all = TRUE)
row.names(ramipril.data.cleaned) <- ramipril.data.cleaned$Geneid
ramipril.data.cleaned$Geneid <- NULL

ramipril_GSE226353_Processed <- ramipril.data.cleaned
new_column_names <- c(rep("0", 3), rep("1", 3))
names(ramipril_GSE226353_Processed) <- new_column_names

# Convert Ensembl IDs to gene symbols using AnnotationDbi
library(AnnotationDbi)
library(org.Mm.eg.db)

ensembl_ids <- row.names(ramipril_GSE226353_Processed)
gene_symbols <- mapIds(org.Mm.eg.db, keys = ensembl_ids, keytype = "ENSEMBL", column = "SYMBOL")

# Create a data frame with Ensembl IDs and gene symbols
gene_mapping <- data.frame(ensembl_id = ensembl_ids, gene_symbol = gene_symbols)

# Remove rows with NA gene symbols
gene_mapping <- gene_mapping[!is.na(gene_mapping$gene_symbol), ]

# Keep only unique gene symbols
gene_mapping <- gene_mapping[!duplicated(gene_mapping$gene_symbol), ]

# Merge the gene mapping with the processed data frame
ramipril_GSE226353_Processed <- merge(gene_mapping, ramipril_GSE226353_Processed, by.x = "ensembl_id", by.y = 0, all.y = TRUE)

# Remove rows with duplicate gene symbols
ramipril_GSE226353_Processed <- ramipril_GSE226353_Processed[!duplicated(ramipril_GSE226353_Processed$gene_symbol), ]
ramipril_GSE226353_Processed <- ramipril_GSE226353_Processed[!is.na(ramipril_GSE226353_Processed$gene_symbol), ]

# Set gene symbols as row names
row.names(ramipril_GSE226353_Processed) <- ramipril_GSE226353_Processed$gene_symbol

# Remove the Ensembl ID and gene symbol columns
ramipril_GSE226353_Processed$ensembl_id <- NULL
ramipril_GSE226353_Processed$gene_symbol <- NULL

# Write on File
# Define the file path
file_path <- "Processed/ramipril_GSE226353_Processed.txt"

# Use write.table to save the matrix to a file
write.table(ramipril_GSE226353_Processed,
            file = file_path, sep = "\t", quote = FALSE,
            row.names = TRUE)


# Canagliflozin -----------------------------------------------------------

canagliflozin.file.data <-  "Canagliflozin_SGLT2_GSE106156/GSE106156_MeanSignals.csv"

canagliflozin.data.raw <- read_csv(canagliflozin.file.data)

# Cleaning Data

canagliflozin.data <- as.data.frame(canagliflozin.data.raw)
canagliflozin.data.cleaned <- na.omit(canagliflozin.data, subset = "Genes")

canagliflozin.data.cleaned <- canagliflozin.data.cleaned %>%
  distinct(Genes, .keep_all = TRUE)
row.names(canagliflozin.data.cleaned) <- canagliflozin.data.cleaned$Genes
canagliflozin.data.cleaned$Genes <- NULL
# rami.df <- cbind(canagliflozin.data.cleaned[,c(5:7)], 
#                 canagliflozin.data.cleaned[,c(2:4)])
Canagliflozin_GSE106156_Processed <- canagliflozin.data.cleaned
new_column_names <- c(rep("0", 4), rep("1", 4))
names(Canagliflozin_GSE106156_Processed) <- new_column_names

# Write on File 

# Define the file path
file_path <- "Processed/Canagliflozin_GSE106156_Processed.txt"

# Use write.table to save the matrix to a file
write.table(Canagliflozin_GSE106156_Processed, 
            file = file_path, sep = "\t", quote = FALSE, 
            row.names = TRUE)



# Ramipril, Furosemide, Hydrocortisone, Cortisone, Valsartan, Losa --------

# Ramipril, Furosemide, Hydrocortisone, Cortisone, Valsartan, 
# Losartan, Prednisolone, Prednisone, Spironolactone

library(GEOquery)
toxicity.data.all <- getGEO("GSE59913")

gene_symbols <- toxicity.data.all[["GSE59913-GPL5426_series_matrix.txt.gz"]]@featureData@data[["GENE_SYMBOL"]]
gene_symbols <- gsub("_predicted", "", gene_symbols)

drugs <- c("Saline", "Ramipril", "Furosemide", "Hydrocortisone", 
           "Cortisone", "Valsartan", "Losartan", "Prednisolone", 
           "Prednisone", "Spironolactone")

get_indices <- function(tissue, drug, day) {
  idx <- grep(paste0("^", tissue, ", ", drug, ", ", day, "$"), 
              toxicity.data.all[["GSE59913-GPL5426_series_matrix.txt.gz"]]@phenoData@data[["source_name_ch1"]])
  return(idx)
}


indices <- unlist(lapply(drugs, function(drug) get_indices("Kidney", drug, "5d")))
indices14 <- unlist(lapply(drugs, function(drug) get_indices("Kidney", drug, "14d")))
expression_data <- toxicity.data.all[["GSE59913-GPL5426_series_matrix.txt.gz"]]@assayData[["exprs"]][, indices]

column_names <- toxicity.data.all[["GSE59913-GPL5426_series_matrix.txt.gz"]]@phenoData@data[["source_name_ch1"]][indices]

master_df <- data.frame(expression_data)
colnames(master_df) <- column_names
# rownames(master_df) <- gene_symbols
master_df$genes <- gene_symbols

library(dplyr)
master_df_cleaned <- master_df %>%
  distinct(genes, .keep_all = TRUE)
master_df_cleaned <- as.data.frame(master_df_cleaned)
row.names(master_df_cleaned) <- master_df_cleaned$genes
master_df_cleaned$genes <- NULL
master_df_cleaned <- master_df_cleaned[-c(8),]
master.df <- subset(master_df_cleaned, select = -c(1:15) )
master.df <- na.omit(master.df)
control_indices <- grep("Kidney, Saline, 5d", column_names)
control_indices <- as.integer(c(1:3))


furosemide <- as.data.frame(master.df[,c(1:3,6:8)])
colnames(furosemide) <- c(rep(0, 3), rep(1, 3))

hydrocortisone <- as.data.frame(master.df[,c(1:3,9:11)])
colnames(hydrocortisone) <- c(rep(0, 3), rep(1, 3))

cortisone <- as.data.frame(master.df[,c(1:3,12:14)])
colnames(cortisone) <- c(rep(0, 3), rep(1, 3))

valsartan <- as.data.frame(master.df[,c(1:3,15:17)])
colnames(valsartan) <- c(rep(0, 3), rep(1, 3))

losartan <- as.data.frame(master.df[,c(1:3,18:20)])
colnames(losartan) <- c(rep(0, 3), rep(1, 3))

prednisolone <- as.data.frame(master.df[,c(1:3,21:23)])
colnames(prednisolone) <- c(rep(0, 3), rep(1, 3))

prednisone <- as.data.frame(master.df[,c(1:3,24:26)])
colnames(prednisone) <- c(rep(0, 3), rep(1, 3))

spironolactone <- as.data.frame(master.df[,c(1:3,27:29)])
colnames(spironolactone) <- c(rep(0, 3), rep(1, 3))

# Save data frames as separate files
folder_path <- "/Users/ahmadali/Library/CloudStorage/OneDrive-Personal/BrandeisUniversity/__Research_120__/Data/Processed"
folder_path <- "C:/Users/ahmad/Library/CloudStorage/OneDrive-Personal/BrandeisUniversity/__Research_120__/Data/Processed"

folder_path <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# write.table(ramipril, file.path(folder_path, "Ramipril_GSE59913_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(furosemide, file.path(folder_path, "Processed/Furosemide_GSE59913_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(hydrocortisone, file.path(folder_path, "Processed/Hydrocortisone_GSE59913_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cortisone, file.path(folder_path, "Processed/Cortisone_GSE59913_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(valsartan, file.path(folder_path, "Processed/Valsartan_GSE59913_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(losartan, file.path(folder_path, "Processed/Losartan_GSE59913_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(prednisolone, file.path(folder_path, "Processed/Prednisolone_GSE59913_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(prednisone, file.path(folder_path, "Processed/Prednisone_GSE59913_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(spironolactone, file.path(folder_path, "Processed/Spironolactone_GSE59913_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)


furosemide_df <- as.data.frame(master.df[,c(1:3,6:8)])
colnames(furosemide_df) <- new_column_names
row.names(furosemide_df) <- row.names(master.df)

write.table(furosemide_df, file.path(folder_path, "Furosemide_GSE59913_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)





# Lisinopril --------------------------------------------------------------
# GSE199437_DKD_all_rpkm_counts.csv.gz


lisinopril.file.data <-  "Lisinopril_GSE199437/GSE199437_DKD_all_rpkm_counts.csv.gz"
lisinopril.data.raw <- read_csv(lisinopril.file.data)
lisinopril.data.raw <- as.data.frame(lisinopril.data.raw)

# Extract gene symbols (first column)
gene_symbols <- lisinopril.data.raw[, 1]

# Create lisinopril dataframe
lisinopril_df <- lisinopril.data.raw[, c(1, 12:21, 22:31)]
colnames(lisinopril_df) <- c("gene_symbol", paste0("0", 1:10), paste0("1", 1:10))

lisinopril_df <- lisinopril_df %>%
  distinct(gene_symbol, .keep_all = TRUE)
lisinopril_df <- as.data.frame(lisinopril_df)
row.names(lisinopril_df) <- lisinopril_df$gene_symbol
lisinopril_df$gene_symbol <- NULL

# Create the other drug dataframe
sglt2i_df <- lisinopril.data.raw[, c(1, 12:21, 32:41)]
colnames(sglt2i_df) <- c("gene_symbol", paste0("0", 1:10), paste0("1", 1:10))

sglt2i_df <- sglt2i_df %>%
  distinct(gene_symbol, .keep_all = TRUE)
sglt2i_df <- as.data.frame(sglt2i_df)
row.names(sglt2i_df) <- sglt2i_df$gene_symbol
sglt2i_df$gene_symbol <- NULL
new_column_names <- c(rep("0", 10), rep("1", 10))
names(lisinopril_df) <- new_column_names
names(sglt2i_df) <- new_column_names

write.table(lisinopril_df, file.path(folder_path, "Processed/Lisinopril_GSE199437_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(sglt2i_df, file.path(folder_path, "Processed/JNJ39933673_GSE199437_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)



# Toxicity 2nd Group ------------------------------------------------------

library(GEOquery)
toxicity.data.all <- getGEO("GSE59913")

gene_symbols <- toxicity.data.all[["GSE59913-GPL5425_series_matrix.txt.gz"]]@featureData@data[["GENE_SYMBOL"]]
gene_symbols <- gsub("_predicted", "", gene_symbols)

drugs <- c("Saline", "Lovastatin", "Simvastatin", "Ibuprofen", "Naproxen", 
           "Methotrexate", "Candesartan", "Lisinopril")

get_indices <- function(tissue, drug, day) {
  idx <- grep(paste0("^", tissue, ", ", drug, ", ", day, "$"), 
              toxicity.data.all[["GSE59913-GPL5425_series_matrix.txt.gz"]]@phenoData@data[["source_name_ch1"]])
  return(idx)
}


indices <- unlist(lapply(drugs, function(drug) get_indices("Kidney", drug, "3d")))
# indices14 <- unlist(lapply(drugs, function(drug) get_indices("Kidney", drug, "14d")))
expression_data <- toxicity.data.all[["GSE59913-GPL5425_series_matrix.txt.gz"]]@assayData[["exprs"]][, indices]

column_names <- toxicity.data.all[["GSE59913-GPL5425_series_matrix.txt.gz"]]@phenoData@data[["source_name_ch1"]][indices]

master_df <- data.frame(expression_data)
colnames(master_df) <- column_names
# rownames(master_df) <- gene_symbols
master_df$genes <- gene_symbols

library(dplyr)
master_df_cleaned <- master_df %>%
  distinct(genes, .keep_all = TRUE)
master_df_cleaned <- as.data.frame(master_df_cleaned)
row.names(master_df_cleaned) <- master_df_cleaned$genes
master_df_cleaned$genes <- NULL
master_df_cleaned <- master_df_cleaned[-c(8),]
# master.df <- subset(master_df_cleaned, select = -c(1:15) )
master.df <- na.omit(master_df_cleaned)
control_indices <- grep("Kidney, Saline, 3d", column_names)

Ibuprofen <- as.data.frame(master.df[,c(control_indices,6:11)])
colnames(Ibuprofen) <- c(rep(0, 5), rep(1, 6))

Naproxen <- as.data.frame(master.df[,c(control_indices,12:17)])
colnames(Naproxen) <- c(rep(0, 5), rep(1, 6))

Methotrexate <- as.data.frame(master.df[,c(control_indices,18:23)])
colnames(Methotrexate) <- c(rep(0, 5), rep(1, 6))

Candesartan <- as.data.frame(master.df[,c(control_indices,24:25)])
colnames(Candesartan) <- c(rep(0, 5), rep(1, 2))

Lisinopril <- as.data.frame(master.df[,c(control_indices,26:28)])
colnames(Lisinopril) <- c(rep(0, 5), rep(1, 3))


folder_path <- setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


write.table(Ibuprofen, file.path(folder_path, "Processed/Ibuprofen_GSE59913_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Naproxen, file.path(folder_path, "Processed/Naproxen_GSE59913_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Methotrexate, file.path(folder_path, "Processed/Methotrexate_GSE59913_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Candesartan, file.path(folder_path, "Processed/Candesartan_GSE59913_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(Lisinopril, file.path(folder_path, "Processed/Lisinopril_GSE59913_Processed.txt"), sep = "\t", quote = FALSE, row.names = TRUE)




