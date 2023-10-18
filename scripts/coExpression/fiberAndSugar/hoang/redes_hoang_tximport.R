# Reset R variables
rm(list = ls())

# Read DESeq2 tutorial
vignette("DESeq2")

##### Configure directory #####

# Notebook
HOME_DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/hoang"

# PC CENA
HOME_DIR = "/home/felipevzps/Documentos/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/hoang"
setwd(HOME_DIR)

# List files in home directory
list.files(HOME_DIR)

###############################

##### Import files: samples, tx2gene file, quantification matrix #####

# Read samples file 
samples <- read.table(file.path(HOME_DIR, "metadata_toCollapse.txt"), header = TRUE) #samples.txt

#Set quant.sf files
files <- file.path(HOME_DIR, "smallData", samples$run, "quant.sf")
all(file.exists(files))

# Set tx2gene file (clusters from MMSeqs2)
tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptomeClassificationTable_0.8_smallData.tsv"), header = FALSE, sep = "\t")
tx2gene

# Organize columns for tx2gene format (transcript ID     group)
tx2gene <- tx2gene[, c(3,2)]
tx2gene

library(tximport)

# Import quantification matrix with tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

names(txi)
head(txi$counts)

######################################################################

##### Calculate the Coefficient of Variation (CV) for each gene #####

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

#####################################################################

library(DESeq2)

# Create DESeqDataSet from txi
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ 1)
dds

##### collapseReplicates (tratando como se fossem technical replicates) #####

?collapseReplicates

ddsColl <- collapseReplicates(dds, dds$sample, dds$run)

# examine the colData and column names of the collapsed data
colData(ddsColl)
ddsColl

##### collapseReplicates (tratando como se fossem technical replicates) #####

##### Pre-filtering #####

# In the tutorial they removed rows that have at least 10 reads total

keep <- rowSums(counts(ddsColl)) >= 1
dds <- ddsColl[keep,]
dds
# 921 rows left

# I want to remove rows with more than 80% of zeros

# Calcular a proporção de zeros em cada linha
zero_prop <- rowSums(counts(dds) == 0) / ncol(counts(dds))

# Defina um limite de 80% para zeros
threshold <- 0.80

# Selecione linhas com menos de 80% de zeros
keep <- zero_prop <= threshold
dds <- dds[keep,]
dds
# 8 rows left

#########################

# Run variance stabilizing transformation on the counts
object <- vst(dds)
object

?vst

# Nao consegui aplicar o vst hoje... (18/10)
# Using shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(ddsColl, normalized=FALSE) + 1),
                           colData=colData(ddsColl))
se

# Plot PCA

# the call to DESeqTransform() is needed to trigger our plotPCA method.
plotPCA( DESeqTransform( se ),intgroup="run" )

pcaData <- plotPCA( DESeqTransform( se ),intgroup="run", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, shape=run, color = run)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

##### Effects of transformations on the variance #####

# shifted logarithm transformation - log2(n + 1)
ntd <- normTransform(dds)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("vsn")

library("vsn")

meanSdPlot(assay(ntd))

# variance stabilizing transformation

meanSdPlot(assay(object))

##### Effects of transformations on the variance #####

### plot PCA

?plotPCA 
plotPCA(dds, intgroup=c("sample", "runsCollapsed"))

###
# Trying PCA
# Code from biostars.org/p/9560363

# Calculate the variance for each gene
# The assay function is used to extract the matrix of normalized values.
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
