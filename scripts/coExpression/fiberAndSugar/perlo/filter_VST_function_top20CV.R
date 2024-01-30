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
HOME_DIR = "/Storage/data1/felipe.peres/Sugarcane_ncRNA/9_Fiber_and_Sugar/co-expression/Perlo/code/coding"
setwd(HOME_DIR)

# *** 1 - Import samples, VST expression matrix, tx2gene, cv ***

# Read samples file 
samples <- read.table(file.path(HOME_DIR, 'infos_perlo_metadata.tsv'), header = TRUE, sep = '\t')

vst_matrix <- read.table(file.path(HOME_DIR, "Perlo2022_counts_filters_VST.txt"))

# Set tx2gene file (clusters from OrthoFinder and MMSeqs2)
tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv", header = FALSE, sep = "\t"))
#tx2gene <- read.table(file.path(HOME_DIR, "panTranscriptome_panRNAomeClassificationTable_hyphen_Class_smallData.tsv"), header = FALSE, sep = "\t")
print("tx2gene file (clusters from MMSeqs2 + OrthoFinder)")
tx2gene

cv <- read.table(file.path(HOME_DIR, "Perlo2022_counts_filters_VST_CV.txt"))

# *** 3 - Filter matrix by function (coding or non-coding) ***

# Filter VST matrix for protein-coding genes
coding_genes <- intersect(rownames(vst_matrix), tx2gene$V2[tx2gene$V4 %in% c("protein-coding", "protein and non-coding")])
vst_matrix_coding <- vst_matrix[coding_genes, ]
#write.table(vst_matrix_coding, file = file.path(HOME_DIR, "Perlo2022_counts_filters_VST_coding.txt"), sep = "\t", quote = FALSE)

# Filter VST matrix for non-coding genes
noncoding_genes <- intersect(rownames(vst_matrix), tx2gene$V2[tx2gene$V4 == "non-coding"])
vst_matrix_noncoding <- vst_matrix[noncoding_genes, ]
#write.table(vst_matrix_noncoding, file = file.path(HOME_DIR, "Perlo2022_counts_filters_VST_noncoding.txt"), sep = "\t", quote = FALSE)

# *** 4 - Filter matrix by top 20% CV *** 

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
write.table(vst_matrix_coding_top, file = file.path(HOME_DIR, "Perlo2022_counts_filters_VST_coding_top20CV.txt"), sep = "\t", quote = FALSE)
write.table(vst_matrix_noncoding_top, file = file.path(HOME_DIR, "Perlo2022_counts_filters_VST_noncoding_top20CV.txt"), sep = "\t", quote = FALSE)

# *** 5 - Plot PCA (coding) ***

# Calculate PCA with all genes using prcomp
pca_result <- prcomp(t(vst_matrix_coding_top), scale. = TRUE)

# Get PC scores
pca_scores <- as.data.frame(pca_result$x)

# *** Adicionar informações de genótipo e internode type ao DataFrame ***

# *** Adicione uma coluna ao seu DataFrame de amostras indicando o internode_type ***
samples$internode_type <- sub("^.*?_(\\S+)_\\d+-weeks.*$", "\\1", samples$Run)
#samples$internode_type

# *** Adicione uma coluna ao seu DataFrame de amostras indicando a replicate ***
#samples$replicate <- sub("^.*?_([^_]+)_\\S+$", "\\1", samples$Run) # Capturing just "5", "8" and "Ex-5"
samples$replicate <- sub("^.*?_([^_]+_\\d+-weeks)_\\S+$", "\\1", samples$Run)
#samples$replicate

# Adicione uma coluna ao seu DataFrame de amostras indicando o nome do genótipo
samples$genotype <- sub("^.*_(\\S+)$", "\\1", samples$Run)
#samples$genotype

print("dds_vst genotypes")
samples$genotype

print("dds_vst replicates")
samples$replicate

print("dds_vst internode type")
samples$internode_type

print("dds_vst samples")
samples$Run

unique_genotypes <- unique(samples$genotype)
unique_genotypes_fixed <- gsub("[-–]", ".", unique_genotypes)
pca_scores$genotype <- rep(unique_genotypes_fixed, each = nrow(pca_scores) / length(unique_genotypes_fixed))

unique_replicate <- unique(samples$replicate)
unique_replicate_fixed <- gsub("[-–]", ".", unique_replicate)
pca_scores$replicate <- rep(unique_replicate_fixed, each = nrow(pca_scores) / length(unique_replicate_fixed))

unique_internode_type <- unique(samples$internode_type)
unique_internode_type_fixed <- gsub("[-–]", ".", unique_internode_type)
pca_scores$internode_type <- rep(unique_internode_type_fixed, each = nrow(pca_scores) / length(unique_internode_type_fixed))
##
pca_scores$genotype <- samples$genotype
pca_scores$replicate <- samples$replicate
pca_scores$internode_type <- samples$internode_type

# *** Plotar PCA usando ggplot2 ***
percentVar <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = pca_scores$replicate, label = pca_scores$genotype)) +
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
  stat_ellipse(geom = "polygon", level=0.95, alpha=0.1, aes(fill = pca_scores$replicate), color=NA, show.legend = FALSE) + # add ellipse with 95% confidence intervals
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
print("saving PCA as: Perlo2022_VST_PCA_withTissues_coding_top20CV.png")
ggsave("Perlo2022_VST_PCA_withTissues_coding_top20CV.png", pca_plot, bg = "white")

# *** 5 - Plot PCA (non-coding) ***

# Calculate PCA with all genes using prcomp
pca_result <- prcomp(t(vst_matrix_noncoding_top), scale. = TRUE)

# Get PC scores
pca_scores <- as.data.frame(pca_result$x)

# *** Adicionar informações de genótipo e internode type ao DataFrame ***

# *** Adicione uma coluna ao seu DataFrame de amostras indicando o internode_type ***
samples$internode_type <- sub("^.*?_(\\S+)_\\d+-weeks.*$", "\\1", samples$Run)
#samples$internode_type

# *** Adicione uma coluna ao seu DataFrame de amostras indicando a replicate ***
#samples$replicate <- sub("^.*?_([^_]+)_\\S+$", "\\1", samples$Run) # Capturing just "5", "8" and "Ex-5"
samples$replicate <- sub("^.*?_([^_]+_\\d+-weeks)_\\S+$", "\\1", samples$Run)
#samples$replicate

# Adicione uma coluna ao seu DataFrame de amostras indicando o nome do genótipo
samples$genotype <- sub("^.*_(\\S+)$", "\\1", samples$Run)
#samples$genotype

print("dds_vst genotypes")
samples$genotype

print("dds_vst replicates")
samples$replicate

print("dds_vst internode type")
samples$internode_type

print("dds_vst samples")
samples$Run

unique_genotypes <- unique(samples$genotype)
unique_genotypes_fixed <- gsub("[-–]", ".", unique_genotypes)
pca_scores$genotype <- rep(unique_genotypes_fixed, each = nrow(pca_scores) / length(unique_genotypes_fixed))

unique_replicate <- unique(samples$replicate)
unique_replicate_fixed <- gsub("[-–]", ".", unique_replicate)
pca_scores$replicate <- rep(unique_replicate_fixed, each = nrow(pca_scores) / length(unique_replicate_fixed))

unique_internode_type <- unique(samples$internode_type)
unique_internode_type_fixed <- gsub("[-–]", ".", unique_internode_type)
pca_scores$internode_type <- rep(unique_internode_type_fixed, each = nrow(pca_scores) / length(unique_internode_type_fixed))
##
pca_scores$genotype <- samples$genotype
pca_scores$replicate <- samples$replicate
pca_scores$internode_type <- samples$internode_type

# *** Plotar PCA usando ggplot2 ***
percentVar <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = pca_scores$replicate, label = pca_scores$genotype)) +
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
  stat_ellipse(geom = "polygon", level=0.95, alpha=0.1, aes(fill = pca_scores$replicate), color=NA, show.legend = FALSE) + # add ellipse with 95% confidence intervals
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
print("saving PCA as: Perlo2022_VST_PCA_withTissues_noncoding_top20CV.png")
ggsave("Perlo2022_VST_PCA_withTissues_noncoding_top20CV.png", pca_plot, bg = "white")