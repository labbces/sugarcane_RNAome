library(ggplot2)
library(ggrepel)
library(viridisLite)
library(viridis)
library(tximport)
library(DESeq2)

# Reset R variables
##rm(list = ls())

# Read DESeq2 tutorial
#vignette("DESeq2")

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
print("List files in home directory")
list.files(HOME_DIR)

###############################

##### Import files: samples, tx2gene file, quantification matrix #####

# Read samples file 
samples <- read.table(file.path(HOME_DIR, 'infos_hoang_metadata.tsv'), header = TRUE, skip = 1, sep = '\t')

#Set quant.sf files
files <- file.path(HOME_DIR, "../data", samples$Accession, "quant.sf")
print("All file exists")
all(file.exists(files))

# Set tx2gene file (clusters from MMSeqs2)
tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptomeClassificationTable_0.8.tsv"), header = FALSE, sep = "\t")
print("tx2gene file (clusters from MMSeqs2)")
tx2gene

# Organize columns for tx2gene format (transcript ID     group)
tx2gene <- tx2gene[, c(3,2)]
print("New tx2gene -> transcript ID    group")
tx2gene

# Import quantification matrix with tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

print("names txi")
names(txi)

print("head txi")
head(txi$counts)

######################################################################

##### Calculate the Coefficient of Variation (CV) for each gene #####

cv <- apply(txi$counts, 1, function(x) sd(x) / mean(x) * 100)

# Add the CV as a new column to txi object
txi$cv <- cv
names(txi)
print(txi$cv)

# Saving txi file with CV
write.table(txi, file = "QuantificationMatrix_CoefficientVariation.tsv", sep = "\t", row.names = FALSE)

# Open a PNG device for saving the plot
png(filename = "QuantificationMatrix_CoefficientVariation.png", width = 800, height = 600)

# Plot a histogram of the Coefficient of Variation (CV)
hist(txi$cv, breaks = 100, main = "Coefficient of Variation Distribution",
     xlab = "Coefficient of Variation (%)", ylab = "Frequency")

# Close the PNG device to save the plot
dev.off()

#####################################################################

# Create DESeqDataSet from txi
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ 1)

print("dds object")
dds

##### collapseReplicates (tratando como se fossem technical replicates) #####

##?collapseReplicates

ddsColl <- collapseReplicates(dds, dds$Run, dds$Accession)

# examine the colData and column names of the collapsed data
print("columns from collapsed dds (ddsColl)")
colData(ddsColl)

print("complete ddsColl object")
ddsColl

##### collapseReplicates (tratando como se fossem technical replicates) #####

# Calcular a proporção de zeros em cada linha
zero_prop <- rowSums(counts(dds) == 0) / ncol(counts(dds))

# Defina um limite de 80% para zeros
threshold <- 0.80

# Selecione linhas com menos de 80% de zeros
keep <- zero_prop <= threshold
keep_ddsColl <- ddsColl[keep,]

print("keep_ddsColl object")
keep_ddsColl

#########################

##### Calculate the Coefficient of Variation (CV) for each gene #####

print('calculating cv after zeros removal ...')
cv_after_zeros_removal <- apply(counts(keep_ddsColl), 1, function(x) sd(x) / mean(x) * 100)

# Add the CV as a new column to ddsColl object
colData(keep_ddsColl)$cv <- cv_after_zeros_removal

# Open a PNG device for saving the plot
png(filename = "QuantificationMatrix_CoefficientVariation_afterZerosRemoval.png", width = 800, height = 600)

# Plot a histogram of the Coefficient of Variation (CV)
hist(cv_after_zeros_removal, breaks = 50, main = "Coefficient of Variation Distribution after Zeros Removal",
     xlab = "Coefficient of Variation (%)", ylab = "Frequency")

# Close the PNG device to save the plot
dev.off()

##### Remove genes with less than 140% cv

dds_after_cv_filter <- ddsColl[rownames(keep_ddsColl)[cv_after_zeros_removal >= 140],]
dds_after_cv_filter

colData(dds_after_cv_filter)
counts(dds_after_cv_filter)

##### calculate cv again after 140% cv removal
##counts(dds_after_cv_filter)
##dim(counts(dds_after_cv_filter))

##cv_after_cv_removal <- apply(counts(dds_after_cv_filter), 1, function(x) sd(x) / mean(x) * 100)
##cv_after_cv_removal

# Open a PNG device for saving the plot
##png(filename = "QuantificationMatrix_CoefficientVariation_afterZerosRemoval_afterLowCVRemoval.png", width = 800, height = 600)

# Plot a histogram of the Coefficient of Variation (CV)
##hist(cv_after_cv_removal, breaks = 50, main = "Coefficient of Variation Distribution after Zeros Removal and Low CV Removal",
##     xlab = "Coefficient of Variation (%)", ylab = "Frequency")

# Close the PNG device to save the plot
##dev.off()

###### calculate cv again after 140% cv removal

##### Run variance stabilizing transformation on the counts ####

# Adicione um pseudocount de 1 a todas as contagens
pseudocount <- 1
dds_counts <- counts(dds_after_cv_filter)
dds_counts_pseudo <- dds_counts + pseudocount

# Crie um novo objeto DESeqDataSet com as contagens ajustadas
dds_pseudo <- DESeqDataSetFromMatrix(countData = dds_counts_pseudo,
                                     colData = colData(dds_after_cv_filter),
                                     design = ~ 1)

# Execute a transformação VST diretamente no objeto DESeqDataSet
dds_vst <- varianceStabilizingTransformation(dds_pseudo)
#-- note: fitType='parametric', but the dispersion trend was not well captured by the
#function: y = a/x + b, and a local regression fit was automatically substituted.
#specify fitType='local' or 'mean' to avoid this message next time.

# Insira o pseudocount na tabela de transformação
dds_vst$pseudocount <- pseudocount

print("dds_vst samples")
dds_vst$Accession

##### Plotar PCA com VST #####

pca_plot_withoutTissues <- plotPCA( DESeqTransform( dds_vst ),intgroup="Accession" )
ggsave("plot_pca_vst_withoutTissues.png", pca_plot_withoutTissues, bg = "white")

##### Plot diferenciando tecidos top e bottom #####

colors <- viridis::viridis(40) #40 cores

# Adicione uma coluna ao seu DataFrame de amostras indicando se é top ou bottom
dds_vst$internode_type <- sub(".*_(top|bottom)-internode$", "\\1", dds_vst$Run)

# Adicione uma coluna ao seu DataFrame de amostras indicando o nome do genótipo
dds_vst$genotype <- sub("^(.*?)_.*", "\\1", dds_vst$Run)
dds_vst$genotype

print("dds_vst internode types")
dds_vst$internode_type

print("dds_vst samples")
dds_vst$Run

# Extrair matriz de counts ajustadas apos VST
counts_matrix_vst <- assay(dds_vst)

# As colunas sao as samples colapsadas e as linhas sao os grupos de non-coding
# Salvar a matriz em um arquivo CSV

print("saving vst without 80 zeros matrix")
#write.table(counts_matrix_vst, file = "Hoang2017_tpm_vst_without80zeros.txt", sep = "\t", quote = FALSE)

# Plot PCA usando ggplot2 para personalização adicional
pca_data <- plotPCA(DESeqTransform(dds_vst), intgroup = "internode_type", returnData = TRUE)

# Calculando a variância explicada por PC1 e PC2
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

# Saving PCA
print("saving PCA as: plot_pca_vst_withTissues.png")
ggsave("plot_pca_vst_withTissues.png", pca_plot, bg = "white")

################################################
##### Plotar PCA sem VST #####

# Using shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(dds_pseudo, normalized=FALSE) + 1),
                           colData=colData(dds_pseudo))
#se

# Plot PCA

# the call to DESeqTransform() is needed to trigger our plotPCA method.
pca_plot_se <- plotPCA( DESeqTransform( se ),intgroup="sample" )
ggsave("plot_pca_SE_withoutTissues.png", pca_plot_se, bg = "white")

##########################
