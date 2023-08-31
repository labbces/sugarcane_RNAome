#read raw matrix
HOME_DIR = "/home/felipe/Documents/SugarcaneCoExpression/hoang"
setwd(HOME_DIR)

#rm(list = ls())

# expression matrix in TPM
raw <- read.table("head1000000_hoang_merged_quant_counts.txt", head = T, row.names = 1)

#convertendo valores em inteiros
raw <- round(raw)

row_data <- read.table("metadata.txt")
colnames(raw) <- row_data$V2

#> install.packages(c("HiClimR", "wesanderson", "WGCNA"))

library(dplyr)
library(tidyr)

## make a raw PCA for quality control, in this part we don't want to aggregate replicates
library(ggplot2)
library(ggrepel)
library(DESeq2) 
library(WGCNA) #instalar
library(tibble)


raw.good <- raw[goodGenes(t(raw)),]
raw.vst <- vst(as.matrix(raw.good))

raw.pca <- prcomp(x = t(raw.vst), scale = TRUE)
var_explained <- raw.pca$sdev^2/sum(raw.pca$sdev^2)

library(viridisLite)
library(viridis)
colors <- viridis::viridis(20) #20 cores

pca_plot <- raw.pca$x %>% 
  as.data.frame %>%
  rownames_to_column("group") %>%
  separate(group,c("genotype", "tissue"), sep = "_") %>%
  unite(reps, c("genotype", "tissue"), remove = F)  %>% 
  ggplot(aes(x=PC1,y=PC2, label = tissue, fill = genotype )) + 
  geom_text(check_overlap = TRUE, vjust = 1.5) +
  geom_point(aes(color= genotype, shape = tissue ), size = 4 ) + 
  scale_color_manual(values = colors) +
  theme_bw(base_size=20) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top",
        text = element_text(family = "Times New Roman", size=20), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(colour="black", angle = 8, vjust = 0.7, hjust=0.5),
        axis.text.y=element_text(colour="black"),
        axis.line = element_line(colour = "black"))

# Salve o PCA como um arquivo PNG
ggsave("pca_plot.png", pca_plot, width = 20, height = 16, units = "in", dpi = 300)
