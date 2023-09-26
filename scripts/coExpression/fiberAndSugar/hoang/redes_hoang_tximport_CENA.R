# Reset R variables
#rm(list = ls())

# Configure directory
HOME_DIR = "/Storage/data1/felipe.peres/Sugarcane_ncRNA/9_Fiber_and_Sugar/co-expression/Hoang/code"
setwd(HOME_DIR)

# List files in home directory
list.files(HOME_DIR)

# Read samples file 
samples <- read.table(file.path(HOME_DIR, "samples.txt"), header = TRUE)
samples

# Set quant.sf files 
files <- file.path(HOME_DIR, "../data", samples$run, "quant.sf")
all(file.exists(files))

# Set tx2gene file (clusters from MMSeqs2)
tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptomeClassificationTable_0.8.tsv"), header = FALSE, sep = "\t")
head(tx2gene)

# Organize columns for tx2gene format (transcript ID     group)
tx2gene <- tx2gene[, c(3,2)]
head(tx2gene)

library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)
head(txi$counts)

# Calculate the Coefficient of Variation (CV) for each gene
cv <- apply(txi$counts, 1, function(x) sd(x) / mean(x) * 100)

# Add the CV as a new column to txi object
txi$cv <- cv
names(txi)
head(txi$cv)

# Saving txi file with CV
write.table(txi, file = "txi_with_cv.tsv", sep = "\t", row.names = FALSE)

# Open a PNG device for saving the plot
png(filename = "cv_histogram.png", width = 800, height = 600)

# Plot a histogram of the Coefficient of Variation (CV)
hist(txi$cv, breaks = 20, main = "Coefficient of Variation Distribution",
     xlab = "Coefficient of Variation (%)", ylab = "Frequency")

# Close the PNG device to save the plot
dev.off()
