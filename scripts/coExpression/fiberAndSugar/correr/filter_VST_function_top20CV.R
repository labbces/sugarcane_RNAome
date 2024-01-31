library(ggplot2)
library(ggrepel)

# *** Pipeline ***
# 1 - Import samples, VST expression matrix, tx2gene, cv
# 2 - Filter matrix by function (coding or non-coding)
# 3 - Filter matrix by top 20% CV
# 4 - Plot PCA (coding and non-coding)

# *** Reset R variables ***

#rm(list = ls())

# *** Configure directory ***

# My laptop
#HOME_DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/correr"

# PC CENA
#HOME_DIR = "/home/felipevzps/Documentos/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/correr"
#setwd(HOME_DIR)

# Cluster
HOME_DIR = "/Storage/data1/felipe.peres/Sugarcane_ncRNA/9_Fiber_and_Sugar/co-expression/Correr/code/updated_filters/CNC"
setwd(HOME_DIR)

# *** 1 - Import samples, VST expression matrix, tx2gene, cv ***

# Read samples file 
samples <- read.table(file.path(HOME_DIR, 'infos_correr_metadata.tsv'), header = TRUE, sep = '\t')

vst_matrix <- read.table(file.path(HOME_DIR, "Correr2020_counts_filters_VST.txt"))

# Set tx2gene file (clusters from OrthoFinder and MMSeqs2)
tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv"), header = FALSE, sep = "\t")
#tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptome_panRNAomeClassificationTable_hyphen_Class_smallData.tsv"), header = FALSE, sep = "\t")
print("tx2gene file (clusters from MMSeqs2 + OrthoFinder)")
tx2gene

cv <- read.table(file.path(HOME_DIR, "Correr2020_counts_filters_VST_CV.txt"))

# *** 2 - Filter matrix by function (coding or non-coding) ***

# Filter VST matrix for protein-coding genes
coding_genes <- intersect(rownames(vst_matrix), tx2gene$V2[tx2gene$V4 %in% c("protein-coding", "protein and non-coding")])
vst_matrix_coding <- vst_matrix[coding_genes, ]
#write.table(vst_matrix_coding, file = file.path(HOME_DIR, "Correr2020_counts_filters_VST_coding.txt"), sep = "\t", quote = FALSE)

# Filter VST matrix for non-coding genes
noncoding_genes <- intersect(rownames(vst_matrix), tx2gene$V2[tx2gene$V4 == "non-coding"])
vst_matrix_noncoding <- vst_matrix[noncoding_genes, ]
#write.table(vst_matrix_noncoding, file = file.path(HOME_DIR, "Correr2020_counts_filters_VST_noncoding.txt"), sep = "\t", quote = FALSE)

# *** 3 - Filter matrix by top 20% CV *** 

genes_cv <- cv$V1
cv_values <- cv$V2

# Keep 20% protein and non-coding
num_genes_to_keep_coding <- round(0.2 * nrow(vst_matrix_coding))
num_genes_to_keep_noncoding <- round(0.2 * nrow(vst_matrix_noncoding))

# Get top genes (top CV for coding and non-coding)
top_genes_coding <- head(genes_cv[genes_cv %in% rownames(vst_matrix_coding)][order(cv_values[genes_cv %in% rownames(vst_matrix_coding)], decreasing = TRUE)], num_genes_to_keep_coding)
vst_matrix_coding_top <- vst_matrix_coding[top_genes_coding, ]

top_genes_noncoding <- head(genes_cv[genes_cv %in% rownames(vst_matrix_noncoding)][order(cv_values[genes_cv %in% rownames(vst_matrix_noncoding)], decreasing = TRUE)], num_genes_to_keep_noncoding)
vst_matrix_noncoding_top <- vst_matrix_noncoding[top_genes_noncoding, ]

# Save top 20% VST expression matrix
write.table(vst_matrix_coding_top, file = file.path(HOME_DIR, "Correr2020_counts_filters_VST_coding_top20CV.txt"), sep = "\t", quote = FALSE)
write.table(vst_matrix_noncoding_top, file = file.path(HOME_DIR, "Correr2020_counts_filters_VST_noncoding_top20CV.txt"), sep = "\t", quote = FALSE)

# *** 4 - Plot PCA (coding) ***

# Calculate PCA with all genes using prcomp
pca_result <- prcomp(t(vst_matrix_coding_top), scale. = TRUE)

# Get PC scores
pca_scores <- as.data.frame(pca_result$x)

# Add biomass column
samples$biomass_group <- sub(".*_(high|low)", "\\1", samples$Run)
samples$biomass_group

# Add genotype column
samples$genotype <- sub("^(.*?)_.*", "\\1", samples$Run)
samples$genotype

print("ddsColl_top_20_percent genotypes")
samples$genotype

print("ddsColl_top_20_percent internode types")
samples$biomass_group

print("ddsColl_top_20_percent samples")
samples$Run

unique_genotypes <- unique(samples$genotype)
unique_genotypes_fixed <- gsub("[-–]", ".", unique_genotypes)
pca_scores$genotype <- rep(unique_genotypes_fixed, each = nrow(pca_scores) / length(unique_genotypes_fixed))

unique_biomass_groups <- unique(samples$biomass_group)
unique_biomass_groups_fixed <- gsub("[-–]", ".", unique_biomass_groups)
pca_scores$biomass_group <- rep(unique_biomass_groups_fixed, each = nrow(pca_scores) / length(unique_biomass_groups_fixed))

# Plot PCA using ggplot2
percentVar <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = pca_scores$biomass_group, label = pca_scores$genotype)) +
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
  stat_ellipse(geom = "polygon", level=0.95, alpha=0.1, aes(fill = pca_scores$biomass_group), color=NA, show.legend = FALSE) + # add ellipse with 95% confidence intervals
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
print("saving PCA as: Correr2020_VST_PCA_withTissues_coding_top20CV.png")
ggsave("Correr2020_VST_PCA_withTissues_coding_top20CV.png", pca_plot, bg = "white")

# *** 5 - Plot PCA (non-coding) ***

# Calculate PCA with all genes using prcomp
pca_result <- prcomp(t(vst_matrix_noncoding_top), scale. = TRUE)

# Get PC scores
pca_scores <- as.data.frame(pca_result$x)

# Add biomass column
samples$biomass_group <- sub(".*_(high|low)", "\\1", samples$Run)
samples$biomass_group

# Add genotype column
samples$genotype <- sub("^(.*?)_.*", "\\1", samples$Run)
samples$genotype

print("ddsColl_top_20_percent genotypes")
samples$genotype

print("ddsColl_top_20_percent internode types")
samples$biomass_group

print("ddsColl_top_20_percent samples")
samples$Run

unique_genotypes <- unique(samples$genotype)
unique_genotypes_fixed <- gsub("[-–]", ".", unique_genotypes)
pca_scores$genotype <- rep(unique_genotypes_fixed, each = nrow(pca_scores) / length(unique_genotypes_fixed))

unique_biomass_groups <- unique(samples$biomass_group)
unique_biomass_groups_fixed <- gsub("[-–]", ".", unique_biomass_groups)
pca_scores$biomass_group <- rep(unique_biomass_groups_fixed, each = nrow(pca_scores) / length(unique_biomass_groups_fixed))

# Plot PCA using ggplot2
percentVar <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = pca_scores$biomass_group, label = pca_scores$genotype)) +
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
  stat_ellipse(geom = "polygon", level=0.95, alpha=0.1, aes(fill = pca_scores$biomass_group), color=NA, show.legend = FALSE) + # add ellipse with 95% confidence intervals
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
print("saving PCA as: Correr2020_VST_PCA_withTissues_noncoding_top20CV.png")
ggsave("Correr2020_VST_PCA_withTissues_noncoding_top20CV.png", pca_plot, bg = "white")
