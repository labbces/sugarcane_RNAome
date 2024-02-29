library(ggplot2)
library(ggrepel)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript filter_VST_function_topCV.R cv_threshold")
}

# *** Pipeline ***
# 1 - Import samples, VST expression matrix, tx2gene, cv
# 2 - Filter matrix by function (CNC, coding and non-coding)
# 3 - Filter matrix by top 20% CV
# 4 - Plot PCA (CNC, coding and non-coding)

# *** Reset R variables ***

#rm(list = ls())

# *** Configure directory ***

# My laptop
#HOME_DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/perlo"

# PC CENA
#HOME_DIR = "/home/felipevzps/Documentos/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/perlo"

# Cluster
HOME_DIR = "/Storage/data1/felipe.peres/Sugarcane_ncRNA/9_Fiber_and_Sugar/co-expression/Perlo/code/updated_filters/CNC"
setwd(HOME_DIR)

# *** 1 - Import samples, VST expression matrix, tx2gene, cv ***

# Read samples file 
samples <- read.table(file.path(HOME_DIR, 'infos_perlo_metadata.tsv'), header = TRUE, sep = '\t')

#vst_matrix <- read.table(file.path(HOME_DIR, "10k_Perlo2022_counts_filters_VST.txt"), check.names = FALSE)
vst_matrix <- read.table(file.path(HOME_DIR, "Perlo2022_counts_filters_VST.txt"), check.names = FALSE)

# Set tx2gene file (clusters from OrthoFinder and MMSeqs2)
tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv"), header = FALSE, sep = "\t")
#tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptome_panRNAomeClassificationTable_hyphen_Class_smallData.tsv"), header = FALSE, sep = "\t")
print("tx2gene file (clusters from MMSeqs2 + OrthoFinder)")
tx2gene

#cv <- read.table(file.path(HOME_DIR, "10k_Perlo2022_counts_filters_VST_CV.txt"))
cv <- read.table(file.path(HOME_DIR, "Perlo2022_counts_filters_VST_CV.txt"))

# *** 2 - Filter matrix by function (CNC, coding and non-coding) ***

# Filter VST matrix for CNC genes
CNC_genes <- intersect(rownames(vst_matrix), tx2gene$V2[tx2gene$V4 %in% c("protein-coding", "protein and non-coding", "non-coding")])
vst_matrix_CNC <- vst_matrix[CNC_genes, ]
#write.table(vst_matrix_CNC, file = file.path(HOME_DIR, "Perlo2022_counts_filters_VST_CNC.txt"), sep = "\t", quote = FALSE)

# Filter VST matrix for protein-coding genes
coding_genes <- intersect(rownames(vst_matrix), tx2gene$V2[tx2gene$V4 %in% c("protein-coding", "protein and non-coding")])
vst_matrix_coding <- vst_matrix[coding_genes, ]
#write.table(vst_matrix_coding, file = file.path(HOME_DIR, "Perlo2022_counts_filters_VST_coding.txt"), sep = "\t", quote = FALSE)

# Filter VST matrix for non-coding genes
noncoding_genes <- intersect(rownames(vst_matrix), tx2gene$V2[tx2gene$V4 == "non-coding"])
vst_matrix_noncoding <- vst_matrix[noncoding_genes, ]
#write.table(vst_matrix_noncoding, file = file.path(HOME_DIR, "Perlo2022_counts_filters_VST_noncoding.txt"), sep = "\t", quote = FALSE)

# *** 3 - Filter matrix by CV > 2.0 *** 

genes_cv <- cv$V1
cv_values <- cv$V2
threshold_cv <- as.numeric(args[1])
#threshold_cv <- 1

# Filter genes with CV > 2.0 for CNC, coding and non-coding
top_genes_CNC <- genes_cv[cv_values > threshold_cv & genes_cv %in% rownames(vst_matrix_CNC)]
vst_matrix_CNC_top <- vst_matrix_CNC[top_genes_CNC, ]

top_genes_coding <- genes_cv[cv_values > threshold_cv & genes_cv %in% rownames(vst_matrix_coding)]
vst_matrix_coding_top <- vst_matrix_coding[top_genes_coding, ]

top_genes_noncoding <- genes_cv[cv_values > threshold_cv & genes_cv %in% rownames(vst_matrix_noncoding)]
vst_matrix_noncoding_top <- vst_matrix_noncoding[top_genes_noncoding, ]

# Save filtered VST expression matrix
output_matrix_CNC <- paste0("Perlo2022_counts_filters_VST_CNC_CV_above_", threshold_cv, ".txt")
output_matrix_coding <- paste0("Perlo2022_counts_filters_VST_coding_CV_above_", threshold_cv, ".txt")
output_matrix_noncoding <- paste0("Perlo2022_counts_filters_VST_noncoding_CV_above_", threshold_cv, ".txt")

write.table(vst_matrix_CNC_top, file = file.path(HOME_DIR, output_matrix_CNC), sep = "\t", quote = FALSE)
write.table(vst_matrix_coding_top, file = file.path(HOME_DIR, output_matrix_coding), sep = "\t", quote = FALSE)
write.table(vst_matrix_noncoding_top, file = file.path(HOME_DIR, output_matrix_noncoding), sep = "\t", quote = FALSE)

# *** 4 - Plot PCA (CNC) ***

# Calculate PCA with all genes using prcomp
pca_result <- prcomp(t(vst_matrix_CNC_top), scale. = TRUE)

# Get PC scores
pca_scores <- as.data.frame(pca_result$x)

pca_scores$genotype <- sub(".*_([^_]+)$", "\\1", rownames(pca_scores))
pca_scores$internode_type <- sub("^.*(Internode(?:_\\d+)?(?:_Ex-\\d+)?)_.*", "\\1", rownames(pca_scores))
pca_scores$replicate <- sub(".*_(\\d+-weeks)_.*", "\\1", rownames(pca_scores))

pca_scores$stage <- paste(pca_scores$internode_type, pca_scores$replicate, sep = " ")

# *** Plotar PCA usando ggplot2 ***
percentVar <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = pca_scores$stage, label = pca_scores$genotype)) +
  geom_point(size = 2) +
  geom_text_repel(
    box.padding = 0.1, point.padding = 0.1,
    segment.color = "black", segment.size = 0.1, segment.alpha = 0.5, max.overlaps = 15,
    size = 3, color = "black" # Definir a cor do texto do rótulo como preto
  ) +
  labs(title = "PCA - Perlo2022 Contrasting Genotypes in Fiber and Sugar (CNC)",
       x = paste0("PC1 ", "(", percentVar[1], "%)"),
       y = paste0("PC2 ", "(", percentVar[2], "%)"),
       color = "Stage") +
  stat_ellipse(geom = "polygon", level=0.95, alpha=0.1, aes(fill = pca_scores$stage), color=NA, show.legend = FALSE) + # add ellipse with 95% confidence intervals
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
output_filename <- paste0("Perlo2022_VST_PCA_withTissues_CNC_", threshold_cv, ".png")
print("Saving CNC PCA Tissues")
ggsave(output_filename, pca_plot, bg = "white")

# *** 4 - Plot PCA (coding) ***

# Calculate PCA with all genes using prcomp
pca_result <- prcomp(t(vst_matrix_coding_top), scale. = TRUE)

# Get PC scores
pca_scores <- as.data.frame(pca_result$x)

pca_scores$genotype <- sub(".*_([^_]+)$", "\\1", rownames(pca_scores))
pca_scores$internode_type <- sub("^.*(Internode(?:_\\d+)?(?:_Ex-\\d+)?)_.*", "\\1", rownames(pca_scores))
pca_scores$replicate <- sub(".*_(\\d+-weeks)_.*", "\\1", rownames(pca_scores))

pca_scores$stage <- paste(pca_scores$internode_type, pca_scores$replicate, sep = " ")

# *** Plotar PCA usando ggplot2 ***
percentVar <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = pca_scores$stage, label = pca_scores$genotype)) +
  geom_point(size = 2) +
  geom_text_repel(
    box.padding = 0.1, point.padding = 0.1,
    segment.color = "black", segment.size = 0.1, segment.alpha = 0.5, max.overlaps = 15,
    size = 3, color = "black" # Definir a cor do texto do rótulo como preto
  ) +
  labs(title = "PCA - Perlo2022 Contrasting Genotypes in Fiber and Sugar (protein-coding)",
       x = paste0("PC1 ", "(", percentVar[1], "%)"),
       y = paste0("PC2 ", "(", percentVar[2], "%)"),
       color = "Stage") +
  stat_ellipse(geom = "polygon", level=0.95, alpha=0.1, aes(fill = pca_scores$stage), color=NA, show.legend = FALSE) + # add ellipse with 95% confidence intervals
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
output_filename <- paste0("Perlo2022_VST_PCA_withTissues_coding_", threshold_cv, ".png")
print("Saving coding PCA Tissues")
ggsave(output_filename, pca_plot, bg = "white")

# *** 4 - Plot PCA (non-coding) ***

# Calculate PCA with all genes using prcomp
pca_result <- prcomp(t(vst_matrix_noncoding_top), scale. = TRUE)

# Get PC scores
pca_scores <- as.data.frame(pca_result$x)

pca_scores$genotype <- sub(".*_([^_]+)$", "\\1", rownames(pca_scores))
pca_scores$internode_type <- sub("^.*(Internode(?:_\\d+)?(?:_Ex-\\d+)?)_.*", "\\1", rownames(pca_scores))
pca_scores$replicate <- sub(".*_(\\d+-weeks)_.*", "\\1", rownames(pca_scores))

pca_scores$stage <- paste(pca_scores$internode_type, pca_scores$replicate, sep = " ")

# *** Plotar PCA usando ggplot2 ***
percentVar <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = pca_scores$stage, label = pca_scores$genotype)) +
  geom_point(size = 2) +
  geom_text_repel(
    box.padding = 0.1, point.padding = 0.1,
    segment.color = "black", segment.size = 0.1, segment.alpha = 0.5, max.overlaps = 15,
    size = 3, color = "black" # Definir a cor do texto do rótulo como preto
  ) +
  labs(title = "PCA - Perlo2022 Contrasting Genotypes in Fiber and Sugar (non-coding)",
       x = paste0("PC1 ", "(", percentVar[1], "%)"),
       y = paste0("PC2 ", "(", percentVar[2], "%)"),
       color = "Stage") +
  stat_ellipse(geom = "polygon", level=0.95, alpha=0.1, aes(fill = pca_scores$stage), color=NA, show.legend = FALSE) + # add ellipse with 95% confidence intervals
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
output_filename <- paste0("Perlo2022_VST_PCA_withTissues_noncoding_", threshold_cv, ".png")
print("Saving non-coding PCA Tissues")
ggsave(output_filename, pca_plot, bg = "white")
