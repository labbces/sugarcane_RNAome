library(tximport)
library(readr)

#rm(list = ls())

DIR="/Storage/data1/felipe.peres/Sugarcane_ncRNA/11_lncRNASequenceConservation/GenomicReads/PLOT"

setwd(DIR)

samples <- read.table(file.path(DIR, "samples.txt"), header = TRUE, sep = '\t')
files <- file.path(DIR, "../QUANT", samples$samples, "quant.sf")
all(file.exists(files))

# read data with tximport (using EffectiveLength instead of Length)
txi <- tximport(files, type = "none", txIn = TRUE, txOut = TRUE, 
                txIdCol = "Name", abundanceCol = "TPM", 
                countsCol = "NumReads", lengthCol = "EffectiveLength",
                importer = function(x) readr::read_tsv(x))

sample_names <- samples$samples                   # extract sample names 
colnames(txi$counts) <- sample_names              # add sample names to columns

counts <- txi$counts                              # extract counts
genes <- rownames(counts)                         # extract gene names  
effective_length <- rowMeans(txi$length)          # extract effective length means

# combining quant.sf files
combined_matrix <- data.frame(
  Name = genes,
  S._barberi = rowSums(counts[, grepl("S._barberi", colnames(counts))]),
  S._officinarum = rowSums(counts[, grepl("S._officinarum", colnames(counts))]),
  S._spontaneum = rowSums(counts[, grepl("S._spontaneum", colnames(counts))]),
  S._bicolor = rowSums(counts[, grepl("S._bicolor", colnames(counts))]),
  EffectiveLength = effective_length
)

write.table(combined_matrix, file = "combined_matrix.sf", sep = "\t", quote = FALSE, row.names = FALSE)