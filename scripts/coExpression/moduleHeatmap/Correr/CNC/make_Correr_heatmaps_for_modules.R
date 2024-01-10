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

# *** Read filtered VST matrix ***
vst_path <- "Correr2020_counts_filters_VST_top20CV.txt"
# Define first column as index
vst <- read.table(vst_path, header = TRUE, row.names = 1) #encoding = "UTF-8", check.names = FALSE

# *** Read metadata ***
metadata_path <-"infos_correr_metadata.tsv"
metadata <- read.table(metadata_path, sep = "\t", header = T)

# *** Make vectors for each interesting module ***
for (i in Nmods$'Mod No'){
  assign(paste0("Module", i), modules %>% filter(module_No == i))
}

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

# *** Make vectors for each module ***
for (i in Nmods$'Mod No'){
  assign(paste0("module", i, "_dat"), vst[rownames(eval(as.name(paste0("Module",i)))),])
}

# *** Calculate mean expression for each module ***
for (i in Nmods$'Mod No'){
  assign(paste0("dat", i),as.data.frame(t(colMeans(eval(as.name(paste0("module",i, "_dat"))))))) 
  labelling <-get(paste0("dat", i, sep = ""))
  labels <- paste0("Module ", i, sep = "")
  rownames(labelling) <- labels
  assign(paste0("dat", i, sep= ""), labelling)
}

# *** Create a list to aggregate all modules *** 
dat_list <- list()

for (i in Nmods$'Mod No'){
  dat_list[[i]] <- eval(as.name(paste0("dat",i)))
}

# *** Create a DataFrame from all modules list ***
df <- do.call("rbind", dat_list)

# *** Reorder the rows of 'anot' based on the order of 'Genotype' ***
merged_df <- merge(anot, data.frame(Genotype = colnames(df)), by = "Genotype", all.x = TRUE)

anot <- merged_df[order(match(colnames(df), merged_df$Genotype)), ]

# Update 'anot' rownames -> Genotype + Condition in names
rownames(anot) <- colnames(df)

# *** Plot heatmap with mean values for column (i.e. mean for each module)
#png("meanOf30ModulesHeatmap.png", res = 300, width = 10*800, height = 10*2850)
png("meanOf30ModulesHeatmap.png", res = 300, width = 800, height = 2850)

pheatmap(df,
         main ="Contrasting Genotypes (CNC Modules)" ,
         scale = "row",
         annotation_col = anot,
         show_rownames = T,
         col = colorRampPalette(brewer.pal(11, "RdYlGn"))(299),
         cluster_cols = T,
         cluster_rows = T,
         cellwidth = NA,
         cellheight = 8,
         angle_col = 0)
dev.off()
