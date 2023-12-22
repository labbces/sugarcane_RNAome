library(ggplot2)
library(ggrepel)
library(viridisLite)
library(viridis)
library(tximport)
library(DESeq2)

# *** Pipeline
# *** 1 - Run Variance Stabilizing Transformation on the counts
# *** 2 - Plot Coefficient of Variation (CV) for each gene
# *** 3 - Remove degraded samples
# *** 4 - Remove samples with more than 80% zeros
# *** 5 - Plot CV for good samples
# *** 6 - Filter samples by top X% CV
# *** 7 - Plot CV of filtered samples (genes with most variance - top X%)
# *** 8 - Plot PCA

# *** Reset R variables ***
#rm(list = ls())

# *** Configure directory ***

# *** My laptop ***
HOME_DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/hoang"

# *** PC CENA ***
#HOME_DIR = "/home/felipevzps/Documentos/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/hoang"
#setwd(HOME_DIR)

# *** Configure directory ***
HOME_DIR = "/Storage/data1/felipe.peres/Sugarcane_ncRNA/9_Fiber_and_Sugar/co-expression/Hoang/code"
setwd(HOME_DIR)

# *** List files in home directory ***
print("List files in home directory")
list.files(HOME_DIR)

# *** Import files: samples, tx2gene file, quantification matrix ***

# *** Read samples file *** 
samples <- read.table(file.path(HOME_DIR, 'infos_hoang_metadata.tsv'), header = TRUE, skip = 1, sep = '\t')

# *** Set quant.sf files ***
files <- file.path(HOME_DIR, "../data", samples$Accession, "quant.sf")
files <- file.path(HOME_DIR, "smallData", samples$Accession, "quant.sf")
print("All file exists")
all(file.exists(files))

# *** Set tx2gene file (clusters from MMSeqs2) ***
tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptomeClassificationTable_0.8.tsv"), header = FALSE, sep = "\t")
tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptomeClassificationTable_0.8_smallData.tsv"), header = FALSE, sep = "\t")
tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptome_panRNAomeClassificationTable_smallData.tsv"), header = FALSE, sep = "\t")

print("tx2gene file (clusters from MMSeqs2)")
tx2gene

# *** Organize columns for tx2gene format (transcript ID     group) ***
tx2gene <- tx2gene[, c(3,2)]
print("New tx2gene -> transcript ID    group")
tx2gene

# *** Import quantification matrix with tximport ***
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

print("names txi")
names(txi)

print("head txi")
head(txi$counts)

# *** Create DESeqDataSet from txi ***
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ 1)

print("dds object")
dds

# *** collapseReplicates (tratando como se fossem technical replicates) ***

ddsColl <- collapseReplicates(dds, dds$Run, dds$Accession)
ddsColl

# *** Examine the colData and column names of the collapsed data ***
print("columns from collapsed dds (ddsColl)")
colData(ddsColl)

print("complete ddsColl object")
ddsColl

# *** 1 - Run Variance Stabilizing Transformation on the counts

# *** Adicione um pseudocount de 1 a todas as contagens ***
print('adding pseudocounts to dds_pseudo')
pseudocount <- 1
dds_counts <- counts(ddsColl)
dds_counts_pseudo <- dds_counts + pseudocount

# *** Crie um novo objeto DESeqDataSet com as contagens ajustadas ***
dds_pseudo <- DESeqDataSetFromMatrix(countData = dds_counts_pseudo,
                                     colData = colData(ddsColl),
                                     design = ~ 1)

# *** Execute a transformação VST diretamente no objeto DESeqDataSet ***
print('applying vst to dds_pseudo')
dds_vst <- varianceStabilizingTransformation(dds_pseudo)
#-- note: fitType='parametric', but the dispersion trend was not well captured by the
#function: y = a/x + b, and a local regression fit was automatically substituted.
#specify fitType='local' or 'mean' to avoid this message next time.

# *** Insira o pseudocount na tabela de transformação ***
dds_vst$pseudocount <- pseudocount

print("dds_vst samples")
dds_vst$Accession

# *** 2 - Plot Coefficient of Variation (CV) for each gene after VST normalization

# *** Calculate the Coefficient of Variation (CV) 

print('calculating cv after vst transformation ...')
cv_after_vst <- apply(assay(dds_vst), 1, function(x) sd(x) / mean(x) * 100)

# *** Add the CV as a new row to ddsColl object ***
colData(dds_vst)
rowData(dds_vst)$cv <- cv_after_vst

#colData(dds_vst)
#rowData(dds_vst)

# *** Open a PNG device for saving the plot ***
print('saving cv plot after VST normalization to file: QuantificationMatrix_CoefficientVariation_afterVST.png')
png(filename = "QuantificationMatrix_CoefficientVariation_afterVST.png", width = 800, height = 600)

# *** Plot a histogram of the Coefficient of Variation (CV) ***
hist(cv_after_vst, breaks = 50, main = "Coefficient of Variation Distribution after VST normalization",
     xlab = "Coefficient of Variation (%)", ylab = "Frequency")

# *** Close the PNG device to save the plot ***
dev.off()

# *** 3 - Remove degraded samples

print('removing degraded samples')
withoutDegradedSamples_ddsColl <- dds_vst[, dds_vst$X..Trimmed <= 30]
withoutDegradedSamples_ddsColl

# *** Remove genes with more than 80% zeros ***

# *** Calcular a proporção de zeros em cada linha ***
zero_prop <- rowSums(assay(withoutDegradedSamples_ddsColl) == 0) / ncol(assay(withoutDegradedSamples_ddsColl))

# *** Define zeros threshold ***
threshold <- 0.80

# *** 4 - Remove genes with more than 80% zeros
keep <- zero_prop <= threshold
withoutDegradedSamplesAndZeros_ddsColl <- withoutDegradedSamples_ddsColl[keep,]

print("withoutDegradedSamplesAndZeros_ddsColl object")
withoutDegradedSamplesAndZeros_ddsColl

# *** 5 - Plot CV for good samples

# *** Calculate the Coefficient of Variation (CV) after degraded samples and zeros removal ***

print('calculating cv after degraded samples and zeros removal ...')
cv_after_zeros_removal <- apply(assay(withoutDegradedSamplesAndZeros_ddsColl), 1, function(x) sd(x) / mean(x) * 100)

# *** Add the CV as a new row to ddsColl object ***
colData(withoutDegradedSamplesAndZeros_ddsColl)
rowData(withoutDegradedSamplesAndZeros_ddsColl)$cv <- cv_after_zeros_removal

#colData(withoutDegradedSamplesAndZeros_ddsColl)
#rowData(withoutDegradedSamplesAndZeros_ddsColl)

# *** Open a PNG device for saving the plot ***
print('saving cv plot after degraded samples and zeros removal to file: QuantificationMatrix_CoefficientVariation_afterDegradedSamplesAndZerosRemoval.png')
png(filename = "QuantificationMatrix_CoefficientVariation_afterDegradedSamplesAndZerosRemoval.png", width = 800, height = 600)

# *** Plot a histogram of the Coefficient of Variation (CV) ***
hist(cv_after_zeros_removal, breaks = 50, main = "Coefficient of Variation Distribution after Degraded Samples and Zeros Removal",
     xlab = "Coefficient of Variation (%)", ylab = "Frequency")

# *** Close the PNG device to save the plot ***
dev.off()

# *** 6 - Filter genes by top 20% CV

# *** Sort genes indices based on CV ***
#order(rowData(withoutDegradedSamplesAndZeros_ddsColl)$cv, decreasing = FALSE)
sorted_genes <- order(rowData(withoutDegradedSamplesAndZeros_ddsColl)$cv, decreasing = TRUE)

# *** Calculate the index for the top 20% *** 
top_20_percent_index <- round(length(sorted_genes) * 0.01) #1%

# *** Select the top 20% genes ***
top_20_percent_genes <- rownames(withoutDegradedSamplesAndZeros_ddsColl)[sorted_genes[1:top_20_percent_index]]

# *** Filter the ddsColl object to keep only the top 20% genes ***
ddsColl_top_20_percent <- withoutDegradedSamplesAndZeros_ddsColl[top_20_percent_genes, ]
print("Now, ddsColl_top_20_percent contains only the top 20% genes based on CV")
ddsColl_top_20_percent

# *** 7 - Plot CV of filtered samples (genes with most variance - top 20%)

print('calculating cv after keep only top 20% genes based on CV ...')
cv_after_cv_filter <- apply(assay(ddsColl_top_20_percent), 1, function(x) sd(x) / mean(x) * 100)

# *** Add the CV as a new row to ddsColl object ***
#colData(ddsColl_top_20_percent)
rowData(ddsColl_top_20_percent)$cv <- cv_after_cv_filter

#colData(ddsColl_top_20_percent)
#rowData(ddsColl_top_20_percent)

# *** Open a PNG device for saving the plot ***
print('saving cv plot after keeping only the top 20% genes based on CV: QuantificationMatrix_CoefficientVariation_top20CV.png')
png(filename = "QuantificationMatrix_CoefficientVariation_top20CV.png", width = 800, height = 600)

# *** Plot a histogram of the Coefficient of Variation (CV) ***
hist(cv_after_cv_filter, breaks = 50, main = "Coefficient of Variation Distribution - top 20% genes based on CV",
     xlab = "Coefficient of Variation (%)", ylab = "Frequency")

# *** Close the PNG device to save the plot ***
dev.off()

# *** 8 - Plot PCA

# *** Calcular PCA com todos os genes da matriz filtrada usando prcomp ***

pca_result <- prcomp(t(assay(ddsColl_top_20_percent)), scale. = TRUE)

# *** Obter os scores dos componentes principais ***
pca_scores <- as.data.frame(pca_result$x)

# *** Adicionar informações de genótipo e internode type ao DataFrame ***

# *** Adicione uma coluna ao seu DataFrame de amostras indicando se é top ou bottom ***
ddsColl_top_20_percent$internode_type <- sub(".*_(top|bottom)-internode$", "\\1", ddsColl_top_20_percent$Run)

# *** Adicione uma coluna ao seu DataFrame de amostras indicando sugar content ***
ddsColl_top_20_percent$sugar_content <- sub("^.*?_(.*?)_.*", "\\1", ddsColl_top_20_percent$Run)

# *** Adicione uma coluna ao seu DataFrame de amostras indicando o nome do genótipo ***
ddsColl_top_20_percent$genotype <- sub("^(.*?)_.*", "\\1", ddsColl_top_20_percent$Run)

print("ddsColl_top_20_percent genotypes")
ddsColl_top_20_percent$genotype

print("ddsColl_top_20_percent sugar content")
ddsColl_top_20_percent$sugar_content

print("ddsColl_top_20_percent internode types")
ddsColl_top_20_percent$internode_type

print("ddsColl_top_20_percent samples")
ddsColl_top_20_percent$Run

pca_scores$genotype <- ddsColl_top_20_percent$genotype
pca_scores$sugar_content <- ddsColl_top_20_percent$sugar_content
pca_scores$internode_type <- ddsColl_top_20_percent$internode_type

# *** Plotar PCA usando ggplot2 ***
percentVar <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2))

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = ddsColl_top_20_percent$sugar_content, shape = internode_type, label = ddsColl_top_20_percent$genotype)) +
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

# *** Saving PCA ***
#print(pca_plot)
print("saving PCA as: plot_pca_vst_withTissues.png")
ggsave("plot_pca_vst_withTissues.png", pca_plot, bg = "white")