library(ggplot2)
library(ggrepel)
library(viridisLite)
library(viridis)
library(tximport)
library(DESeq2)

# *** Pipeline
# *** 1 - Plot Coefficient of Variation (CV) for each gene
# *** 2 - Remove degraded samples
# *** 3 - Remove samples with more than 80% zeros
# *** 4 - Plot CV for good samples
# *** 5 - Run Variance Stabilizing Transformation on the counts
# *** 6 - Plot CV after VST
# *** 7 - Filter samples by top 20% CV
# *** 8 - Plot CV of filtered samples (genes with most variance - top X%)
# *** 9 - Plot PCA

# *** Reset R variables ***
#rm(list = ls())

# *** Configure directory ***

# *** My laptop ***
#HOME_DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/perlo"

# *** PC CENA ***
#HOME_DIR = "/home/felipevzps/Documentos/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/perlo"
#setwd(HOME_DIR)

# *** Configure directory ***
HOME_DIR = "/Storage/data1/felipe.peres/Sugarcane_ncRNA/9_Fiber_and_Sugar/co-expression/Perlo/code/non-coding"
setwd(HOME_DIR)

# *** List files in home directory ***
print("List files in home directory")
list.files(HOME_DIR)

# *** Import files: samples, tx2gene file, quantification matrix ***

# *** Read samples file *** 
samples <- read.table(file.path(HOME_DIR, 'infos_perlo_metadata.tsv'), header = TRUE, sep = '\t')

# *** Set quant.sf files ***
files <- file.path(HOME_DIR, "../../data", samples$Accession, "quant.sf")
#files <- file.path(HOME_DIR, "smallData", samples$Accession, "quant.sf")
print("All file exists")
all(file.exists(files))

# *** Set tx2gene file (clusters from MMSeqs2) ***
tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv"), header = FALSE, sep = "\t")
#tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptome_panRNAomeClassificationTable_hyphen.tsv"), header = FALSE, sep = "\t")
#tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptomeClassificationTable_0.8_smallData.tsv"), header = FALSE, sep = "\t")
#tx2gene <- read.table(file.path(HOME_DIR, "../../tx2gene/I2.8/smallData/panTranscriptome_panRNAomeClassificationTable_smallData_hyphen.tsv"), header = FALSE, sep = "\t")

print("tx2gene file (clusters from MMSeqs2)")
tx2gene

# *** Organize columns for tx2gene format (transcript ID     group) ***
tx2gene_non_coding <- subset(tx2gene, V4 == "non-coding")
tx2gene_non_coding <- tx2gene_non_coding[, c(3,2)]
print("New tx2gene -> transcript ID    group")
tx2gene_non_coding

# *** Import quantification matrix with tximport ***
txi <- tximport(files, type = "salmon", tx2gene = tx2gene_non_coding)

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

# *** 1 - Plot Coefficient of Variation (CV) for each gene

# *** Calculate the Coefficient of Variation (CV) 

print('calculating raw cv ...')
raw_cv <- apply(assay(ddsColl), 1, function(x) sd(x) / mean(x))

# *** Add the CV as a new row to ddsColl object ***
colData(ddsColl)
rowData(ddsColl)$cv <- raw_cv

#colData(ddsColl)
#rowData(ddsColl)

# *** Open a PNG device for saving the plot ***
print('saving raw cv to file: QuantificationMatrix_RawCoefficientVariation.png')
png(filename = "QuantificationMatrix_RawCoefficientVariation.png", width = 800, height = 600)

# *** Plot a histogram of the Coefficient of Variation (CV) ***
hist(raw_cv, breaks = 50, main = "Raw Coefficient of Variation Distribution",
     xlab = "Coefficient of Variation", ylab = "Frequency")

# *** Plot histogram with log transformation on X axis
#hist(log(raw_cv), breaks = 50, main = "Raw Coefficient of Variation Distribution",
#     xlab = "Log(Coefficient of Variation)", ylab = "Frequency")

# *** Add a text annotation for the count of CVs equal to zero ***
count_zero_cv <- sum(raw_cv == 0)
text(0, 0, sprintf("CVs = 0: %d", count_zero_cv), adj = c(0, 1), col = "red", cex = 1.2)

# *** Close the PNG device to save the plot ***
dev.off()

# *** 2 - Remove degraded samples

print('removing degraded samples')
withoutDegradedSamples_ddsColl <- ddsColl[, ddsColl$X..Trimmed <= 30]
withoutDegradedSamples_ddsColl

# *** 3 - Remove genes with more than 80% zeros ***

# *** Calcular a proporção de zeros em cada linha ***
zero_prop <- rowSums(assay(withoutDegradedSamples_ddsColl) == 0) / ncol(assay(withoutDegradedSamples_ddsColl))

# *** Define zeros threshold ***
threshold <- 0.80

keep <- zero_prop <= threshold
withoutDegradedSamplesAndZeros_ddsColl <- withoutDegradedSamples_ddsColl[keep,]

print("withoutDegradedSamplesAndZeros_ddsColl object")
withoutDegradedSamplesAndZeros_ddsColl

# *** 4 - Plot CV for good samples

# *** Calculate the Coefficient of Variation (CV) after degraded samples and zeros removal ***

print('calculating cv after degraded samples and zeros removal ...')
cv_after_zeros_removal <- apply(assay(withoutDegradedSamplesAndZeros_ddsColl), 1, function(x) sd(x) / mean(x))

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
     xlab = "Coefficient of Variation", ylab = "Frequency")

# *** Plot histogram with log transformation on X axis
#hist(log(cv_after_zeros_removal), breaks = 50, main = "Coefficient of Variation Distribution after Degraded Samples and Zeros Removal",
#     xlab = "Log(Coefficient of Variation)", ylab = "Frequency")

# *** Add a text annotation for the count of CVs equal to zero ***
count_zero_cv <- sum(cv_after_zeros_removal == 0)
text(0, 0, sprintf("CVs = 0: %d", count_zero_cv), adj = c(0, 1), col = "red", cex = 1.2)

# *** Close the PNG device to save the plot ***
dev.off()

# *** 5 - Run Variance Stabilizing Transformation on the counts

# *** Adicione um pseudocount de 1 a todas as contagens ***
print('adding pseudocounts to dds_pseudo')
pseudocount <- 1
dds_counts <- counts(withoutDegradedSamplesAndZeros_ddsColl)
dds_counts_pseudo <- dds_counts + pseudocount

# *** Crie um novo objeto DESeqDataSet com as contagens ajustadas ***
dds_pseudo <- DESeqDataSetFromMatrix(countData = dds_counts_pseudo,
                                     colData = colData(withoutDegradedSamplesAndZeros_ddsColl),
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

# *** 6 - Plot Coefficient of Variation (CV) for each gene after VST normalization

# *** Calculate the Coefficient of Variation (CV) 

print('calculating cv after vst transformation ...')
cv_after_vst <- apply(assay(dds_vst), 1, function(x) sd(x) / mean(x))

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
     xlab = "Coefficient of Variation", ylab = "Frequency")

# *** Plot histogram with log transformation on X axis
#hist(log(cv_after_vst), breaks = 50, main = "Coefficient of Variation Distribution after VST normalization",
#     xlab = "Log(Coefficient of Variation)", ylab = "Frequency")

# *** Add a text annotation for the count of CVs equal to zero ***
count_zero_cv <- sum(cv_after_vst == 0)
text(0, 0, sprintf("CVs = 0: %d", count_zero_cv), adj = c(0, 1), col = "red", cex = 1.2)

# *** Close the PNG device to save the plot ***
dev.off()

# *** 7 - Filter genes by top 20% CV

# *** Sort genes indices based on CV ***
#order(rowData(dds_vst)$cv, decreasing = FALSE)
sorted_genes <- order(rowData(dds_vst)$cv, decreasing = TRUE)

# *** Calculate the index for the top 20% *** 
top_20_percent_index <- round(length(sorted_genes) * 0.20) #20%

# *** Select the top 20% genes ***
top_20_percent_genes <- rownames(dds_vst)[sorted_genes[1:top_20_percent_index]]

# *** Filter the ddsColl object to keep only the top 20% genes ***
ddsColl_top_20_percent <- dds_vst[top_20_percent_genes, ]
print("Now, ddsColl_top_20_percent contains only the top 20% genes based on CV")
ddsColl_top_20_percent

# *** Extrair matriz de counts ajustadas apos VST para os top20% genes ***
counts_matrix_vst_top20 <- assay(ddsColl_top_20_percent)

# *** Salvar a matriz em um arquivo CSV ***
# OBS: As colunas sao as samples colapsadas e as linhas sao os grupos de non-coding

print("saving counts matrix of top20% genes based on CV after VST")
write.table(counts_matrix_vst_top20, file = "Perlo2022_counts_filters_VST_top20CV.txt", sep = "\t", quote = FALSE)

# *** 8 - Plot CV of filtered samples (genes with most variance - top 20%)

print('calculating cv after keep only top 20% genes based on CV ...')
cv_after_cv_filter <- apply(assay(ddsColl_top_20_percent), 1, function(x) sd(x) / mean(x))

print("saving CV of the top20% genes")
write.table(cv_after_cv_filter, file = "Perlo2022_top20CV.txt", sep = "\t", quote = FALSE)

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
     xlab = "Coefficient of Variation", ylab = "Frequency")

# *** Plot histogram with log transformation on X axis
#hist(log(cv_after_cv_filter), breaks = 50, main = "Coefficient of Variation Distribution - top 20% genes based on CV",
#     xlab = "Log(Coefficient of Variation)", ylab = "Frequency")

# *** Add a text annotation for the count of CVs equal to zero ***
count_zero_cv <- sum(cv_after_cv_filter == 0)
text(0, 0, sprintf("CVs = 0: %d", count_zero_cv), adj = c(0, 1), col = "red", cex = 1.2)

# *** Close the PNG device to save the plot ***
dev.off()

# *** 9 - Plot PCA

# *** Calcular PCA com todos os genes da matriz filtrada usando prcomp ***

pca_result <- prcomp(t(assay(ddsColl_top_20_percent)), scale. = TRUE)

# *** Obter os scores dos componentes principais ***
pca_scores <- as.data.frame(pca_result$x)

# *** Adicionar informações de genótipo e internode type ao DataFrame ***

# *** Adicione uma coluna ao seu DataFrame de amostras indicando o internode_type ***
ddsColl_top_20_percent$internode_type <- sub("^.*?_(\\S+)_\\d+-weeks.*$", "\\1", ddsColl_top_20_percent$Run)
#ddsColl_top_20_percent$internode_type

# *** Adicione uma coluna ao seu DataFrame de amostras indicando a replicate ***
#ddsColl_top_20_percent$replicate <- sub("^.*?_([^_]+)_\\S+$", "\\1", ddsColl_top_20_percent$Run) # Capturing just "5", "8" and "Ex-5"
ddsColl_top_20_percent$replicate <- sub("^.*?_([^_]+_\\d+-weeks)_\\S+$", "\\1", ddsColl_top_20_percent$Run)
#ddsColl_top_20_percent$replicate

# Adicione uma coluna ao seu DataFrame de amostras indicando o nome do genótipo
ddsColl_top_20_percent$genotype <- sub("^.*_(\\S+)$", "\\1", ddsColl_top_20_percent$Run)
#ddsColl_top_20_percent$genotype

print("dds_vst genotypes")
ddsColl_top_20_percent$genotype

print("dds_vst replicates")
ddsColl_top_20_percent$replicate

print("dds_vst internode type")
ddsColl_top_20_percent$internode_type

print("dds_vst samples")
ddsColl_top_20_percent$Run

pca_scores$genotype <- ddsColl_top_20_percent$genotype
pca_scores$replicate <- ddsColl_top_20_percent$replicate
pca_scores$internode_type <- ddsColl_top_20_percent$internode_type

# *** Plotar PCA usando ggplot2 ***
percentVar <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = ddsColl_top_20_percent$replicate, label = ddsColl_top_20_percent$genotype)) +
  geom_point(size = 2) +
  geom_text_repel(
    box.padding = 0.1, point.padding = 0.1,
    segment.color = "black", segment.size = 0.1, segment.alpha = 0.5, max.overlaps = 15,
    size = 3, color = "black" # Definir a cor do texto do rótulo como preto
  ) +
  labs(title = "PCA - Perlo2022 Contrasting Genotypes in Fiber and Sugar",
       x = paste0("PC1 ", "(", percentVar[1], "%)"),
       y = paste0("PC2 ", "(", percentVar[2], "%)"),
       color = "Stage") +
  stat_ellipse(geom = "polygon", level=0.95, alpha=0.1, aes(fill = ddsColl_top_20_percent$replicate), color=NA, show.legend = FALSE) + # add ellipse with 95% confidence intervals
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

# *** Saving PCA ***
print(pca_plot)
print("saving PCA as: plot_pca_vst_withTissues.png")
ggsave("plot_pca_vst_withTissues.png", pca_plot, bg = "white")
