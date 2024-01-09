library(tidyverse)
library(DESeq2)
library(pheatmap)
library(viridis)
library(wesanderson)

rm(list=ls())

DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/moduleHeatmap/Correr/CNC"
setwd(DIR)

# Read modules
Nmods <- read.table("all_mods.txt", header = F)
colnames(Nmods) <- "Mod No"

# files
# this repo have only a small subset of vst counts
modules_path <- "Correr2020_counts_filters_VST_top20CV_mcl_I2.0.formated.csv"
vst_path <- "Correr2020_counts_filters_VST_top20CV.txt"
metadata_path <-"infos_correr_metadata.tsv"

modules <- read.table(modules_path, row.names = 1, header = F)
colnames(modules) <- c("module_No")

# read full vst matrix -> definir a primeira coluna como índice
vst <- read.table(vst_path, header = TRUE, row.names = 1)
metadata <- read.table(metadata_path, sep = "\t", header = T)

# Make vectors for each intersting module
for (i in Nmods$'Mod No'){
  assign(paste0("Module", i), modules %>% filter(module_No == i))
}

sample_table <- metadata

sample_table$Condition <- sub(".*_(high|low)", "\\1", metadata$Run)
sample_table$Condition <- as.factor((sample_table$Condition))

sample_table$Genotype <- sub("^(.*?)_.*", "\\1", metadata$Run)
sample_table$Genotype <- as.factor((sample_table$Genotype)) 

sample_table$Group <- as.factor(paste(sample_table$Genotype, '_', sample_table$Condition, sep=''))

#df <- as.data.frame((assay(vst)))
annotation_col <- sample_table

#Remova Duplicatas em annotation_col
unique_annotation_col <- unique(annotation_col[, c("Group", "Genotype", "Condition")])
anot <- select(unique_annotation_col, "Group","Genotype", "Condition")

for (i in Nmods$'Mod No'){
  assign(paste0("module", i, "_dat"), vst[rownames(eval(as.name(paste0("Module",i)))),])
}

for (i in Nmods$'Mod No'){
  assign(paste0("dat", i),as.data.frame(t(colMeans(eval(as.name(paste0("module",i, "_dat"))))))) 
  labelling <-get(paste0("dat", i, sep = ""))
  labels <- paste0("Module ", i, sep = "")
  rownames(labelling) <- labels
  assign(paste0("dat", i, sep= ""), labelling)
}

dat_list <- list()
for (i in Nmods$'Mod No'){
  dat_list[[i]] <- eval(as.name(paste0("dat",i)))
}

df <- do.call("rbind", dat_list)
rownames(anot) <- colnames(df)

# heat map with mean values for column i.e media por modulo
png("all_hp_samples.png", res = 300, width = 10*800, height = 10*2850)

pheatmap(df,
         main ="Nitrogen Modules" ,
         scale = "row",
         annotation_col = anot,
         show_rownames = T,
         col = inferno(299),
         cluster_cols = T,
         cluster_rows = T,
         cellwidth = NA,
         cellheight = 8) 
dev.off()
