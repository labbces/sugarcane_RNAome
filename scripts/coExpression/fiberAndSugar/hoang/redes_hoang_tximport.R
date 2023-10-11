# Reset R variables
rm(list = ls())

# Configure directory
HOME_DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/hoang"

# PC CENA
HOME_DIR = "/home/felipevzps/Documentos/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/hoang"
setwd(HOME_DIR)

# List files in home directory
list.files(HOME_DIR)

# Read samples file 
samples <- read.table(file.path(HOME_DIR, "samples.txt"), header = TRUE)

#Set quant.sf files
files <- file.path(HOME_DIR, "smallData", samples$run, "quant.sf")
all(file.exists(files))

# Set tx2gene file (clusters from MMSeqs2)
tx2gene <- read.table(file.path(HOME_DIR, "tx2gene_smallData.txt"), header = FALSE, sep = "\t")
tx2gene

# Organize columns for tx2gene format (transcript ID     group)
tx2gene <- tx2gene[, c(3,2)]
tx2gene

library(tximport)
# Dont merge top and bottom
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

names(txi)
head(txi$counts)

# Calculate the Coefficient of Variation (CV) for each gene
cv <- apply(txi$counts, 1, function(x) sd(x) / mean(x) * 100)

# Add the CV as a new column to txi object
txi$cv <- cv
names(txi)
print(txi$cv)

# Saving txi file with CV
write.table(txi, file = "txi_with_cv.tsv", sep = "\t", row.names = FALSE)

# Open a PNG device for saving the plot
png(filename = "small_cv_histogram.png", width = 800, height = 600)

# Plot a histogram of the Coefficient of Variation (CV)
hist(txi$cv, breaks = 50, main = "Coefficient of Variation Distribution",
     xlab = "Coefficient of Variation (%)", ylab = "Frequency")

# Close the PNG device to save the plot
dev.off()

# Trying PCA
# Code from biostars.org/p/9560363

library(DESeq2)

# Create DESeqDataSet from txi
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ 1)
dds

# Run variance stabilizing transformation on the counts
object <- vst(dds)
object

# Calculate the variance for each gene
rv <- rowVars(assay(dds))

# Top n genes by variance to keep
ntop <- 500

# Select the ntop genes by variance
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]

# Perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(dds)[select,]))
pca

# Loading for the first two PCs
loadings <- pca$rotation[, seq_len(2)]
loadings
