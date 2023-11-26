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
samples <- read.table(file.path(HOME_DIR, "metadata_toCollapse.txt"), header = TRUE) #samples.txt

#Set quant.sf files
files <- file.path(HOME_DIR, "../data", samples$run, "quant.sf")
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

library(tximport)

# Import quantification matrix with tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

print("names txi")
names(txi)

print("head txi")
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

print("dds object")
dds

##### collapseReplicates (tratando como se fossem technical replicates) #####

##?collapseReplicates

ddsColl <- collapseReplicates(dds, dds$sample, dds$run)

# examine the colData and column names of the collapsed data
print("columns from collapsed dds (ddsColl)")
colData(ddsColl)

print("complete ddsColl object")
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

print("keep_ddsColl object")
keep_ddsColl

#########################

##### Run variance stabilizing transformation on the counts ####

# Adicione um pseudocount de 1 a todas as contagens
pseudocount <- 1
dds_counts <- counts(keep_ddsColl)
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

print("dds_vst samples")
dds_vst$sample

##### Plotar PCA com VST #####
library(ggplot2)

pca_plot <- plotPCA( DESeqTransform( dds_vst ),intgroup="sample" )
ggsave("plot_pca_vst_withoutTissues.png", pca_plot)

##### Plot diferenciando tecidos top e bottom #####

library(viridisLite)
library(viridis)
colors <- viridis::viridis(40) #40 cores

# Adicione uma coluna ao seu DataFrame de amostras indicando se é top ou bottom
dds_vst$internode_type <- sub(".*_(top|bottom)-internode$", "\\1", dds_vst$sample)

print("dds_vst internode types")
dds_vst$internode_type

print("dds_vst samples")
dds_vst$sample

# Extrair matriz de counts ajustadas apos VST
counts_matrix_vst <- assay(dds_vst)

# As colunas sao as samples colapsadas e as linhas sao os grupos de non-coding
# Salvar a matriz em um arquivo CSV

print("saving vst without 80 zeros matrix")
write.csv(counts_matrix_vst, file = "Hoang2017_tpm_vst_without80zeros.txt")

# Plot PCA usando ggplot2 para personalização adicional
pca_data <- plotPCA(DESeqTransform(dds_vst), intgroup = "internode_type", returnData = TRUE)

# Personalize o gráfico
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = dds_vst$sample, shape = internode_type)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot",
       x = paste0("PC1 (", round(100 * pca_data$percentVar[1], 1), "%)"),
       y = paste0("PC2 (", round(100 * pca_data$percentVar[2], 1), "%)")) +
  theme_minimal()

# Adicione uma escala de cores manualmente para corresponder aos shapes
pca_plot <- pca_plot + scale_color_manual(values = colors, name = "Sample")

# salvar o gráfico
ggsave("plot_pca_vst_withTissues.png", pca_plot)

################################################
##### Plotar PCA sem VST #####

# Using shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(ddsColl, normalized=FALSE) + 1),
                           colData=colData(ddsColl))
#se

# Plot PCA

# the call to DESeqTransform() is needed to trigger our plotPCA method.
pca_plot <- plotPCA( DESeqTransform( se ),intgroup="sample" )
ggsave("plot_pca_SE_withoutTissues.png", pca_plot)

##########################
