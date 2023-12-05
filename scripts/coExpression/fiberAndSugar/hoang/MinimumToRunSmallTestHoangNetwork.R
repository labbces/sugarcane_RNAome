library(ggplot2)
library(ggrepel)
library(viridisLite)
library(viridis)
library(tximport)
library(DESeq2)

#Minimum code to run redes_hoang_tximport.R 

# Reset R variables
rm(list = ls())

HOME_DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/hoang"
setwd(HOME_DIR)

samples <- read.table(file.path(HOME_DIR, 'infos_hoang_metadata.tsv'), header = TRUE, skip = 1, sep = '\t')
files <- file.path(HOME_DIR, "smallData", samples$Accession, "quant.sf")

tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptomeClassificationTable_0.8_smallData.tsv"), header = FALSE, sep = "\t")

tx2gene <- tx2gene[, c(3,2)]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

cv <- apply(txi$counts, 1, function(x) sd(x) / mean(x) * 100)
txi$cv <- cv

hist(txi$cv, breaks = 50, main = "Coefficient of Variation Distribution",
     xlab = "Coefficient of Variation (%)", ylab = "Frequency")

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ 1)
ddsColl <- collapseReplicates(dds, dds$Run, dds$Accession)

zero_prop <- rowSums(counts(dds) == 0) / ncol(counts(dds))
threshold <- 0.80

keep <- zero_prop <= threshold
ddsColl <- ddsColl[keep,]

cv_after_zeros_removal <- apply(counts(ddsColl), 1, function(x) sd(x) / mean(x) * 100)
colData(ddsColl)$cv <- cv_after_zeros_removal


hist(cv_after_zeros_removal, breaks = 50, main = "Coefficient of Variation Distribution after Zeros Removal",
     xlab = "Coefficient of Variation (%)", ylab = "Frequency")

ddsColl <- ddsColl[, ddsColl$X..Trimmed <= 30]

dds_after_cv_filter <- ddsColl[rownames(ddsColl)[cv_after_zeros_removal >= 140],]

cv_after_cv_removal <- apply(counts(dds_after_cv_filter), 1, function(x) sd(x) / mean(x) * 100)

hist(cv_after_cv_removal, breaks = 50, main = "Coefficient of Variation Distribution after Zeros Removal and Low CV Removal",
     xlab = "Coefficient of Variation (%)", ylab = "Frequency")

dds_vst <- dds_after_cv_filter

colors <- viridis::viridis(40) #40 cores
dds_vst$internode_type <- sub(".*_(top|bottom)-internode$", "\\1", dds_vst$Run)
dds_vst$genotype <- sub("^(.*?)_.*", "\\1", dds_vst$Run)
pca_data <- plotPCA(DESeqTransform(dds_vst), intgroup = "internode_type", returnData = TRUE)
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

print(pca_plot)

