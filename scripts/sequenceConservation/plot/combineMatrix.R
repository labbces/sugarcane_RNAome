library(tximport)
library(readr)

#rm(list = ls())

#DIR="/home/felipe/Documents/sequenceConservation_test_combine_matrix"
DIR="/Storage/data1/felipe.peres/Sugarcane_ncRNA/11_lncRNASequenceConservation/GenomicReads/PLOT"

setwd(DIR)
#list.files(DIR)

samples <- read.table(file.path(DIR, 'samples.txt'), header = TRUE, sep = '\t')
files <- file.path(DIR, samples$samples, "quant.sf")
all(file.exists(files))

# Read data with tximport (without tx2gene)
txi <- tximport(files, type = "none", txIn = TRUE, txOut = TRUE, 
                txIdCol = "Name", abundanceCol = "TPM", 
                countsCol = "NumReads", lengthCol = "Length",
                importer = function(x) readr::read_tsv(x))

# counts
#head(txi$counts)

# abundance (TPM)
#head(txi$abundance)

# Extract sample names from samples.txt
sample_names <- samples$samples

# Add names to columns
colnames(txi$counts) <- sample_names

# Extract counts and gene names
counts <- txi$counts
#head(counts)

genes <- rownames(counts)

# Combining new DataFrame
combined_matrix <- data.frame(
  Name = genes,
  S._barberi = rowSums(counts[, grepl("S._barberi", colnames(counts))]),
  S._officinarum = rowSums(counts[, grepl("S._officinarum", colnames(counts))]),
  S._spontaneum = rowSums(counts[, grepl("S._spontaneum", colnames(counts))]),
  S._bicolor = rowSums(counts[, grepl("S._bicolor", colnames(counts))])
)

# Saving new DataFrame to combined_matrix.sf
write.table(combined_matrix, file = "combined_matrix.sf", sep = "\t", quote = FALSE, row.names = FALSE)
