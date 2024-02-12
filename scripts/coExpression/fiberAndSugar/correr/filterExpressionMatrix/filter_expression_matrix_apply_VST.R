library(tximport)
library(DESeq2)
library(ggplot2)
library(ggrepel)

# *** Pipeline ***
# 1 - Remove degraded samples (30% + trimmed)
# 2 - Remove genes with 100% zeros
# 3 - Run Variance Stabilizing Transformation on the counts
# 4 - Calculate and plot CV after VST
# 5 - Save expression matrix (filtered + VST)
# 6 - Save CV for each gene (filtered + VST)
# 7 - Plot PCA

# *** Reset R variables ***

#rm(list = ls())

# *** Configure directory ***

# My laptop
#HOME_DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/correr"

# PC CENA
#HOME_DIR = "/home/felipevzps/Documentos/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/correr"

# Cluster
HOME_DIR = "/Storage/data1/felipe.peres/Sugarcane_ncRNA/9_Fiber_and_Sugar/co-expression/Correr/code/updated_filters/CNC/faster"
setwd(HOME_DIR)

# List files in home directory
print("List files in home directory")
list.files(HOME_DIR)

# *** Import files: samples, tx2gene file, expression matrix ***

# Read samples file 
samples <- read.table(file.path(HOME_DIR, 'infos_correr_metadata.tsv'), header = TRUE, sep = '\t')

# Set quant.sf files
files <- file.path(HOME_DIR, "../../../../data", samples$Accession, "quant.sf")
#files <- file.path(HOME_DIR, "smallData", samples$Accession, "quant.sf")

print("All file exists")
all(file.exists(files))

# Set tx2gene file (clusters from OrthoFinder and MMSeqs2)
tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv"), header = FALSE, sep = "\t")
#tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptome_panRNAomeClassificationTable_hyphen_Class_smallData.tsv"), header = FALSE, sep = "\t")
#print("tx2gene file (clusters from MMSeqs2 + OrthoFinder)")
tx2gene

# Organize columns for tx2gene format (transcript ID     gene)
tx2gene <- tx2gene[, c(3,2)]
print("New tx2gene -> transcript ID    gene")
tx2gene

# Import quantification matrix with tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

print("names txi")
names(txi)

print("head txi")
head(txi$counts)

# Create DESeqDataSet from txi
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ 1)

print("dds object")
dds

# collapseReplicates (tratando como se fossem technical replicates)

ddsColl <- collapseReplicates(dds, dds$Run, dds$Accession)
ddsColl

# Examine the colData and column names of the collapsed data
print("columns from collapsed dds (ddsColl)")
colData(ddsColl)

print("complete ddsColl object")
ddsColl

# *** 1 - Remove degraded samples (30% + trimmed) ***

print('removing degraded samples')
withoutDegradedSamples_ddsColl <- ddsColl[, ddsColl$X..Trimmed <= 30]
withoutDegradedSamples_ddsColl

# *** 2 - Remove genes with 100% zeros ***

# Calculate zeros proportion
zero_prop <- rowSums(assay(withoutDegradedSamples_ddsColl) == 0) / ncol(assay(withoutDegradedSamples_ddsColl))

# Define zeros threshold
threshold <- 1

keep <- zero_prop < threshold
withoutDegradedSamplesAndZeros_ddsColl <- withoutDegradedSamples_ddsColl[keep,]

print("withoutDegradedSamplesAndZeros_ddsColl object")
withoutDegradedSamplesAndZeros_ddsColl

# *** 3 - Run Variance Stabilizing Transformation on the counts ***

# Add a pseudocount = 1 to counts
print('adding pseudocounts to dds_pseudo')
pseudocount <- 1
dds_counts <- counts(withoutDegradedSamplesAndZeros_ddsColl)
dds_counts_pseudo <- dds_counts + pseudocount

# Create a new DESeqDataSet object with adjusted counts
dds_pseudo <- DESeqDataSetFromMatrix(countData = dds_counts_pseudo,
                                     colData = colData(withoutDegradedSamplesAndZeros_ddsColl),
                                     design = ~ 1)

# Apply VST directly on DESeqDataSet object
print('applying vst to dds_pseudo')
dds_vst <- varianceStabilizingTransformation(dds_pseudo)
#-- note: fitType='parametric', but the dispersion trend was not well captured by the
#function: y = a/x + b, and a local regression fit was automatically substituted.
#specify fitType='local' or 'mean' to avoid this message next time.

# Insert the pseudocount on transformation table
dds_vst$pseudocount <- pseudocount

print("dds_vst samples")
dds_vst$Accession

# *** 4 - Calculate and plot CV after VST ***

# Calculate the Coefficient of Variation (CV) 
print('calculating cv after vst transformation ...')
cv_after_vst <- apply(assay(dds_vst), 1, function(x) abs(sd(x)) / abs(mean(x)))

# Add the CV as a new row to ddsColl object
colData(dds_vst)
rowData(dds_vst)$cv <- cv_after_vst

# Open a PNG device for saving the plot
print('saving cv plot after VST normalization to file: CNC_ExpressionMatrix_CoefficientVariation_afterVST.png')
png(filename = "CNC_ExpressionMatrix_CoefficientVariation_afterVST.png", width = 800, height = 600)

# Plot a histogram of the Coefficient of Variation (CV)
hist(cv_after_vst, breaks = 50, main = "Coefficient of Variation Distribution after VST normalization (CNC)",
     xlab = "Coefficient of Variation", ylab = "Frequency")

# Plot histogram with log transformation on X axis
#hist(log(cv_after_vst), breaks = 50, main = "Coefficient of Variation Distribution after VST normalization",
#     xlab = "Log(Coefficient of Variation)", ylab = "Frequency")

# Add a text annotation for the count of CVs equal to zero
count_zero_cv <- sum(cv_after_vst == 0)
text(0, 0, sprintf("CVs = 0: %d", count_zero_cv), adj = c(0, 1), col = "red", cex = 1.2)

# Close the PNG device to save the plot
dev.off()

# *** 5 - Save expression matrix (filtered + VST) ***
expression_matrix_vst <- assay(dds_vst)

# Save expression matrix in a csv file
# OBS: Columns are collapsed samples and lines are genes
print("saving expresion matrix after VST")
write.table(expression_matrix_vst, file = "Correr2020_counts_filters_VST.txt", sep = "\t", quote = FALSE)

# *** 6 - Save CV for each gene (filtered + VST)

cv_data <- data.frame(Gene = rownames(dds_vst), CV = cv_after_vst)
cv_data <- cv_data[order(cv_data$CV, decreasing = TRUE), ]

# Save CV for each gene
print("saving CV for each gene")
write.table(cv_data, file = "Correr2020_counts_filters_VST_CV.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# *** 7 - Plot PCA ***

# Calculate PCA with all genes using prcomp
pca_result <- prcomp(t(assay(dds_vst)), scale. = TRUE)

# Get PC scores
pca_scores <- as.data.frame(pca_result$x)

# Add biomass column
dds_vst$biomass_group <- sub(".*_(high|low)", "\\1", dds_vst$Run)
dds_vst$biomass_group

# Add genotype column
dds_vst$genotype <- sub("^(.*?)_.*", "\\1", dds_vst$Run)
dds_vst$genotype

print("ddsColl_top_20_percent genotypes")
dds_vst$genotype

print("ddsColl_top_20_percent internode types")
dds_vst$biomass_group

print("ddsColl_top_20_percent samples")
dds_vst$Run

pca_scores$genotype <- dds_vst$genotype
pca_scores$biomass_group <- dds_vst$biomass_group

# Plot PCA using ggplot2
percentVar <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = dds_vst$biomass_group, label = dds_vst$genotype)) +
  geom_point(size = 2) +
  geom_text_repel(
    box.padding = 0.1, point.padding = 0.1,
    segment.color = "black", segment.size = 0.1, segment.alpha = 0.5, max.overlaps = 15,
    size = 3, color = "black" # Definir a cor do texto do rótulo como preto
  ) +
  labs(title = "PCA - Correr2020 Contrasting Genotypes in Fiber and Sugar",
       x = paste0("PC1 ", "(",percentVar[1], "%)"),
       y = paste0("PC2 ", "(",percentVar[2], "%)"),
       color = "Groups") +
  stat_ellipse(geom = "polygon", level=0.95, alpha=0.1, aes(fill = dds_vst$biomass_group), color=NA, show.legend = FALSE) + # add ellipse with 95% confidence intervals
  theme_classic() +
  theme(
    axis.line = element_blank(),  # Linha dos eixos X e Y
    panel.grid.major = element_line(color = alpha("gray", 0.2)),  # Remover linhas de grade principais
    panel.grid.minor = element_line(color = alpha("gray", 0.2)),  # Remover linhas de grade secundárias
    panel.border = element_rect(color = "transparent", fill = NA),  # Cor da borda do painel
    plot.background = element_rect(fill = "white"),  # Cor do fundo do gráfico
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Adicionar linha pontilhada no eixo x
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")   # Adicionar linha pontilhada no eixo y

# Saving PCA
#print(pca_plot)
print("saving PCA as: Correr2020_CNC_VST_PCA_withTissues.png")
ggsave("Correr2020_CNC_VST_PCA_withTissues.png", pca_plot, bg = "white")
