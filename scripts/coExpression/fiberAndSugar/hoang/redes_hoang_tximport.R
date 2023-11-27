library(ggplot2)
library(ggrepel)
library(viridisLite)
library(viridis)
library(tximport)
library(DESeq2)

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
ddsColl <- ddsColl[keep,]
head(ddsColl)
rownames(ddsColl)
colnames(ddsColl)

gene_de_interesse <- "OG1388542"

counts_gene_interesse <- counts(ddsColl)[gene_de_interesse, ]
print(counts_gene_interesse)

# 8 rows left

#########################

##### Run variance stabilizing transformation on the counts ####

?vst

# Adicione um pseudocount de 1 a todas as contagens
pseudocount <- 1
dds_counts <- counts(ddsColl)
dds_counts_pseudo <- dds_counts + pseudocount

# Crie um novo objeto DESeqDataSet com as contagens ajustadas
dds_pseudo <- DESeqDataSetFromMatrix(countData = dds_counts_pseudo,
                                     colData = colData(ddsColl),
                                     design = ~ 1)

# Execute a transformação VST diretamente no objeto DESeqDataSet
dds_vst <- varianceStabilizingTransformation(dds_pseudo)
#-- note: fitType='parametric', but the dispersion trend was not well captured by the
#function: y = a/x + b, and a local regression fit was automatically substituted.
#specify fitType='local' or 'mean' to avoid this message next time.

# Insira o pseudocount na tabela de transformação
dds_vst$pseudocount <- pseudocount
dds_vst$sample
dds_vst$run
# Extrair a matriz de contagens ajustadas após o VST
counts_matrix_vst <- assay(dds_vst)

# Salvar a matriz 
write.table(counts_matrix_vst, file = "small_counts_matrix_vst.txt", sep = "\t", quote = FALSE)

################################################################

##### Plotando PCA #####

# Just to work with same names in plot (gambiarra temporaria)
dds_vst <- ddsColl

# Plot PCA diferenciando tecidos top e bottom

colors <- viridis::viridis(40) #40 cores

# Adicione uma coluna ao seu DataFrame de amostras indicando se é top ou bottom
dds_vst$internode_type <- sub(".*_(top|bottom)-internode$", "\\1", dds_vst$sample)
dds_vst$internode_type

# Adicione uma coluna ao seu DataFrame de amostras indicando o nome do genótipo
dds_vst$genotype <- sub("^(.*?)_.*", "\\1", dds_vst$sample)
dds_vst$genotype

dds_vst$sample

# Plot PCA usando ggplot2 para personalização adicional
pca_data <- plotPCA(DESeqTransform(dds_vst), intgroup = "internode_type", returnData = TRUE)
print(pca_data)

percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = dds_vst$genotype, shape = internode_type, label = dds_vst$genotype)) +
  geom_point(size = 2) +
  geom_text_repel(
    box.padding = 0.1, point.padding = 0.1,
    segment.color = "black", segment.size = 0.1, segment.alpha = 0.5, max.overlaps = 15,
    size = 2
  ) +
  labs(title = "PCA - Hoang2017 Contrasting Genotypes in Fiber and Sugar",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       color = "Genotype",
       shape = "Internode Type") +
  theme_minimal()

# Exiba o gráfico
print(pca_plot)

ggsave("plot_pca_vst_withTissues.png", pca_plot)

################################################
##### Plotar PCA sem VST #####

# Using shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(ddsColl, normalized=FALSE) + 1),
                           colData=colData(ddsColl))
se

# Plot PCA

# the call to DESeqTransform() is needed to trigger our plotPCA method.
pca_plot <- plotPCA( DESeqTransform( se ),intgroup="sample" )
print(pca_plot)
ggsave("plot_pca_SE_withoutTissues.png", pca_plot)


pcaData <- plotPCA( DESeqTransform( se ),intgroup="run", returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, shape=run, color = run)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

##############################

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

######################################################

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
