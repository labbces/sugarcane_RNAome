# Reset R variables
##rm(list = ls())

# Read DESeq2 tutorial
vignette("DESeq2")

##### Configure directory #####

# Notebook
##HOME_DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/hoang"

# PC CENA
##HOME_DIR = "/home/felipevzps/Documentos/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/hoang"
##setwd(HOME_DIR)

# Configure directory
HOME_DIR = "/Storage/data1/felipe.peres/Sugarcane_ncRNA/9_Fiber_and_Sugar/co-expression/Hoang/code"
setwd(HOME_DIR)

# List files in home directory
list.files(HOME_DIR)

###############################

##### Import files: samples, tx2gene file, quantification matrix #####

# Read samples file 
samples <- read.table(file.path(HOME_DIR, "metadata_toCollapse.txt"), header = TRUE) #samples.txt

#Set quant.sf files
files <- file.path(HOME_DIR, "../data", samples$run, "quant.sf")
all(file.exists(files))

# Set tx2gene file (clusters from MMSeqs2)
tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptomeClassificationTable_0.8.tsv"), header = FALSE, sep = "\t")
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

##cv <- apply(txi$counts, 1, function(x) sd(x) / mean(x) * 100)

# Add the CV as a new column to txi object
##txi$cv <- cv
##names(txi)
##print(txi$cv)

# Saving txi file with CV
##write.table(txi, file = "txi_with_cv.tsv", sep = "\t", row.names = FALSE)

# Open a PNG device for saving the plot
##png(filename = "small_cv_histogram.png", width = 800, height = 600)

# Plot a histogram of the Coefficient of Variation (CV)
##hist(txi$cv, breaks = 50, main = "Coefficient of Variation Distribution",
##     xlab = "Coefficient of Variation (%)", ylab = "Frequency")

# Close the PNG device to save the plot
##dev.off()

#####################################################################

library(DESeq2)

# Create DESeqDataSet from txi
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ 1)
dds

##### collapseReplicates (tratando como se fossem technical replicates) #####

##?collapseReplicates

ddsColl <- collapseReplicates(dds, dds$sample, dds$run)

# examine the colData and column names of the collapsed data
colData(ddsColl)
ddsColl

##### collapseReplicates (tratando como se fossem technical replicates) #####

##### Pre-filtering #####

# In the tutorial they removed rows that have at least 10 reads total

##keep <- rowSums(counts(ddsColl)) >= 1
##dds <- ddsColl[keep,]
##dds
# 921 rows left

# I want to remove rows with more than 80% of zeros

# Calcular a proporção de zeros em cada linha
zero_prop <- rowSums(counts(ddsColl) == 0) / ncol(counts(ddsColl))

# Defina um limite de 80% para zeros
threshold <- 0.80

# Selecione linhas com menos de 80% de zeros
keep <- zero_prop <= threshold
keep_ddsColl <- ddsColl[keep,]
keep_ddsColl
# 8 rows left

#########################
