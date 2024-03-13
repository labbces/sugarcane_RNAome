library(tidyverse)
library(DESeq2)
library(pheatmap)
library(viridis)
library(wesanderson)
library(RColorBrewer)

rm(list=ls())

DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/moduleHeatmap/Correr/CNC"
setwd(DIR)

# *** Read file for modules numbers (1,2,3,4,5 ...) *** 
Nmods <- read.table("all_mods.txt", header = F)
#Nmods <- read.table("one_mods.txt", header = F)
colnames(Nmods) <- "Mod No"

# *** Read formated modules ***
#modules_path <- "Correr2020_counts_filters_VST_top20CV_mcl_I2.0.formated.csv"
# formated cliques
modules_path <- "Correr2020_counts_filters_VST_top20CV_mcl_I2.0.formated_cliques.csv"

# TODO: Check duplicate genes in first column
#modules <- read.table(modules_path, row.names = 1, header = F)
modules <- read.table(modules_path, header = F)
# Verificar e remover linhas duplicadas
modules <- modules[!duplicated(modules$V1), ]
# Atribuir os valores da primeira coluna como nomes de linha
row.names(modules) <- modules$V1
# Remover a primeira coluna dos dados
modules <- modules[, -1]
modules <- modules[, c(2,1)]
# end TODO

#colnames(modules) <- c("module_No")
colnames(modules) <- c("module_No", "classification")

modules$gene <- rownames(modules)

# *** Read filtered VST matrix ***
vst_path <- "Correr2020_counts_filters_VST_CNC_CV_above2.txt"
# Define first column as index
vst <- read.table(vst_path, header = TRUE, row.names = 1, check.names = FALSE) #encoding = "UTF-8", check.names = FALSE

# Extract only the first name before the hyphen in row names
colnames(vst) <- sub("^([^_]+)_.*", "\\1", colnames(vst))
# This regular expression pattern looks for a lowercase letter followed by an uppercase letter and inserts a space between them
colnames(vst) <- sub("([a-z])([A-Z])", "\\1 \\2", colnames(vst))

# *** Read metadata ***
metadata_path <-"infos_correr_metadata.tsv"
metadata <- read.table(metadata_path, sep = "\t", header = T)

# *** Define groups to plot (Groups, Genotypes) ***
sample_table <- metadata

# Groups
sample_table$Groups <- sub(".*_(high|low)", "\\1", metadata$Run)
# Sub hyphen for space
sample_table$Groups <- gsub("-", " ", sample_table$Groups)
sample_table$Groups <- as.factor((sample_table$Groups))

# Genotypes
sample_table$Genotypes <- sub("^(.*?)_.*", "\\1", metadata$Run)
sample_table$Genotypes <- as.factor((sample_table$Genotypes)) 
# This regular expression pattern looks for a lowercase letter followed by an uppercase letter and inserts a space between them
sample_table$Genotypes <- sub("([a-z])([A-Z])", "\\1 \\2", sample_table$Genotypes)

# Group (Genotypes and Groups)
sample_table$Group <- as.factor(paste(sample_table$Genotypes, ' ', sample_table$Groups, sep=''))

annotation_col <- sample_table

# Remove duplicates in annotation_col (metadata)
unique_annotation_col <- distinct(annotation_col, Group, Genotypes, Groups)
anot <- select(unique_annotation_col, Genotypes, Groups)

# red (-) black (0) green (+)
my_palette = colorRampPalette(c("red", "black", "green"))(n=1000)

for (i in Nmods[,1]){
  names <- modules[modules$module_No == i,]
  df <- vst[names$gene,]
  
  # *** Reorder the rows of 'anot' based on the order of 'Genotypes' ***
  merged_df <- merge(anot, data.frame(Genotypes = colnames(df)), by = "Genotypes", all.x = TRUE)
  anot <- merged_df[order(match(colnames(df), merged_df$Genotypes)), ]
  # Force the use of "-" in the rownames instead of "."
  colnames(df) <- gsub("\\.", "-", colnames(df))
  # Update 'anot' rownames -> Genotypes + Groups in names
  rownames(anot) <- colnames(df)
  #colnames(df) <- colnames(vst)
  
  # heat map with mean values for column i.e media por modulo
  png(paste0("module_", i, "_heatmap",".png", sep = ""), res = 300, width = 5*800, height = 5*2850) # Too big
  pheatmap(df,
           main =paste0("Genotypes contrasting in biomass production (CNC Module ",i, ")", sep = "") ,
           scale = "row",
           annotation_col = anot,
           show_rownames = T,
           col = my_palette,
           cluster_cols = T,
           cluster_rows = T,
           cellwidth = NA,
           cellheight = 8,
           angle_col = 45) 
  dev.off()
  print(i)
}
