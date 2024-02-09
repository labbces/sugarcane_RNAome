library(ggplot2)

# *** Reset R variables ***
rm(list = ls())

DIR="/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/hoang/AnalysisCV/top20CV"
DIR="/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/hoang/AnalysisCV/top50CV"

DIR="/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/correr/top20CV"
DIR="/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/correr/AnalysisCV/top50CV"

DIR="/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/perlo/AnalysisCV/top20CV"
DIR="/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/perlo/AnalysisCV/top50CV"

setwd(DIR)

# Leia os coeficientes de variação
genes <- read.table(file.path(DIR,"Hoang2017_top20CV.txt"), header = TRUE)
genes <- read.table(file.path(DIR,"Hoang2017_top50CV.txt"), header = TRUE)

genes <- read.table(file.path(DIR,"Correr2020_top20CV.txt"), header = TRUE)
genes <- read.table(file.path(DIR,"Correr2020_top50CV.txt"), header = TRUE)

genes <- read.table(file.path(DIR,"Perlo2022_top20CV.txt"), header = TRUE)
genes <- read.table(file.path(DIR,"Perlo2022_top50CV.txt"), header = TRUE)

# Leia a matriz de quantificação
quantification_matrix <- read.table(file.path(DIR,'Hoang2017_counts_filters_VST_top20CV.txt'), header = TRUE, row.names = 1, check.names = FALSE)
quantification_matrix <- read.table(file.path(DIR,'Hoang2017_counts_filters_VST_top50CV.txt'), header = TRUE, row.names = 1, check.names = FALSE)

quantification_matrix <- read.table(file.path(DIR,'Correr2020_counts_filters_VST_top20CV.txt'), header = TRUE, row.names = 1, check.names = FALSE)
quantification_matrix <- read.table(file.path(DIR,'Correr2020_counts_filters_VST_top50CV.txt'), header = TRUE, row.names = 1, check.names = FALSE)

quantification_matrix <- read.table(file.path(DIR,'Perlo2022_counts_filters_VST_top20CV.txt'), header = TRUE, row.names = 1, check.names = FALSE)
quantification_matrix <- read.table(file.path(DIR,'Perlo2022_counts_filters_VST_top50CV.txt'), header = TRUE, row.names = 1, check.names = FALSE)

# Divida os genes em quartis
genes$Quartil <- cut(genes$CV, breaks = quantile(genes$CV, probs = seq(0, 1, 0.25), na.rm = TRUE), include.lowest = TRUE, labels = FALSE)

# Analisar os quartis
quantile(genes$CV, probs = seq(0, 1, 0.25), na.rm = TRUE)

#min CV     Q1         Q2          Q3         max CV     
#0%         25%        50%         75%        100% 
#0.5606058  0.5878932  0.6233223   0.6774185  1.5504023
# Ou seja, 75% dos CV estão abaixo de 0.6774185

quartis <- quantile(genes$CV, probs = seq(0, 1, 0.25), na.rm = TRUE)
percent_quartis <- format(quartis, digits = 2)  # Mantendo os valores originais
names(percent_quartis) <- c("Min CV", "25%", "50%", "75%", "Max CV")
print(percent_quartis)

# Define o número de pares de genes desejados por quartil
n_pairs_per_quartil <- 4

# Gere pares de genes aleatórios para cada quartil
random_gene_pairs <- lapply(split(genes, genes$Quartil), function(subset) {
  sample(subset$Gene, size = 2 * n_pairs_per_quartil, replace = TRUE)
})

percent_range <- c("0-25%", "25-50%", "50-75%", "75-100%")

# Extrair as condições das colunas
conditions <- sub(".*_", "", colnames(quantification_matrix))

# Ordenar as colunas de acordo com as condições
ordered_columns <- colnames(quantification_matrix)[order(conditions)]

# Aplicar a nova ordem às colunas
gene_pair_quantification_ordered <- quantification_matrix[, ordered_columns]

for (i in seq_along(random_gene_pairs)) {
  matrix_gene_pairs <- matrix(random_gene_pairs[[i]], ncol = 2, byrow = TRUE)
  
  plots <- list()
  
  for (j in seq_len(n_pairs_per_quartil)) {
    gene_pair <- matrix_gene_pairs[j, ]
    gene_pair_quantification <- gene_pair_quantification_ordered[gene_pair, ]
    
    # Obtenha as cores dinamicamente
    gene_labels <- rownames(gene_pair_quantification)
    
    # Crie uma tabela de dados onde cada linha representa uma combinação de gene e condição
    plot_data <- data.frame(
      Condition = rep(colnames(gene_pair_quantification), each = 1),
      Counts_Gene1 = as.vector(t(gene_pair_quantification[1, ])),
      Counts_Gene2 = as.vector(t(gene_pair_quantification[2, ])),
      Gene_Color1 = gene_labels[1],
      Gene_Color2 = gene_labels[2]
    )
    
    plot_data$Condition <- factor(plot_data$Condition, levels = unique(plot_data$Condition))
    
    # Encontrar os índices onde a condição muda
    mudanca_condicao <- which(diff(grepl("low", plot_data$Condition)) != 0) + 0.5
    
    plot <- ggplot(plot_data, aes(x = plot_data$Condition)) +
      geom_point(aes(y = Counts_Gene1, color = Gene_Color1), position = position_jitter(width = 0.2), shape = 1) +
      geom_point(aes(y = Counts_Gene2, color = Gene_Color2), position = position_jitter(width = 0.2), shape = 2) +
      geom_vline(xintercept = mudanca_condicao, linetype="dashed", color = "red") +  # Adiciona a linha vertical
      
      labs(title = paste("Q", i, "(", percent_range[i], ")", '- Gene Pair', j, "(top 20% CV)"),
           x = 'Samples',
           y = 'Expression Level (VST)',
           color = "Genes") +
      theme_classic() +
      theme(
        axis.line = element_line(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "transparent", fill = NA),
        plot.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5)
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plots[[length(plots) + 1]] <- plot
  }
  
  # Salvar os gráficos em um único arquivo por quartil
  ggsave(paste0("gene_pairs_quartil_", i, ".png"), do.call(gridExtra::arrangeGrob, plots), width = 15, height = 10)
}
