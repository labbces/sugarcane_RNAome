library(tidyverse)
library(DESeq2)
library(pheatmap)
library(viridis)
library(wesanderson)
library(RColorBrewer)

rm(list=ls())

DIR = "/home/felipe/Documents/exploring_Jorge_heatmaps/Correr/CNC/heatmap_per_module"
setwd(DIR)

# *** Read file for modules numbers (1,2,3,4,5 ...) *** 
Nmods <- read.table("all_mods.txt", header = F)
colnames(Nmods) <- "Mod No"

# *** Read formated modules ***
modules_path <- "Correr2020_counts_filters_VST_top20CV_mcl_I2.0.formated.csv"
modules <- read.table(modules_path, row.names = 1, header = F)
colnames(modules) <- c("module_No")
modules$gene <- rownames(modules)

# *** Read filtered VST matrix ***
vst_path <- "Correr2020_counts_filters_VST_top20CV.txt"
# Define first column as index
vst <- read.table(vst_path, header = TRUE, row.names = 1) #encoding = "UTF-8", check.names = FALSE

# *** Read metadata ***
metadata_path <-"infos_correr_metadata.tsv"
metadata <- read.table(metadata_path, sep = "\t", header = T)

# *** Define groups to plot (Condition, Genotype) ***
sample_table <- metadata

# Condition
sample_table$Condition <- sub(".*_(high|low)", "\\1", metadata$Run)
# Sub hyphen for dot
sample_table$Condition <- gsub("-", ".", sample_table$Condition)
sample_table$Condition <- as.factor((sample_table$Condition))

# Genotype
sample_table$Genotype <- sub("^(.*?)_.*", "\\1", metadata$Run)
sample_table$Genotype <- as.factor((sample_table$Genotype)) 

# Group (Genotype and Condition)
sample_table$Group <- as.factor(paste(sample_table$Genotype, '_', sample_table$Condition, sep=''))

annotation_col <- sample_table

# Remove duplicates in annotation_col (metadata)
unique_annotation_col <- distinct(annotation_col, Group, Genotype, Condition)
anot <- select(unique_annotation_col, Genotype, Condition)

for (i in Nmods[,1]){
  names <- modules[modules$module_No == i,]
  df <- vst[names$gene,]
  
  # *** Reorder the rows of 'anot' based on the order of 'Genotype' ***
  merged_df <- merge(anot, data.frame(Genotype = colnames(df)), by = "Genotype", all.x = TRUE)
  anot <- merged_df[order(match(colnames(df), merged_df$Genotype)), ]
  # Update 'anot' rownames -> Genotype + Condition in names
  rownames(anot) <- colnames(df)
  colnames(df) <- colnames(vst)
  
  # heat map with mean values for column i.e media por modulo
  png(paste0("hp_samples_module_", i, ".png", sep = ""), res = 300, width = 5.5*800, height = 5*2850)
  pheatmap(df,
           main =paste0("Contrasting Genotypes (CNC Module ",i, ")", sep = "") ,
           scale = "row",
           annotation_col = anot,
           show_rownames = T,
           col = colorRampPalette(brewer.pal(11, "RdYlGn"))(299),
           cluster_cols = T,
           cluster_rows = T,
           cellwidth = NA,
           cellheight = 8,
           angle_col = 45) 
  dev.off()
  print(i)
}
