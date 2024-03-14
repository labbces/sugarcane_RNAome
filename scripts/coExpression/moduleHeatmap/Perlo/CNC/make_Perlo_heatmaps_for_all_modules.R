library(tidyverse)
library(DESeq2)
library(pheatmap)
library(viridis)
library(wesanderson)
library(RColorBrewer)

rm(list=ls())

DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/moduleHeatmap/Perlo/CNC"
setwd(DIR)

# *** Read file for modules numbers (1,2,3,4,5 ...) *** 
Nmods <- read.table("all_mods.txt", header = F)
colnames(Nmods) <- "Mod No"

# *** Read formated modules ***
#modules_path <- "Correr2020_counts_filters_VST_top20CV_mcl_I2.0.formated.csv"
# formated cliques
modules_path <- "Perlo2022_counts_filters_VST_topCV_mcl_formated_cliques.csv"

# TODO: Check duplicate genes in first column
#modules <- read.table(modules_path, row.names = 1, header = F) # Only works with unique gene names in first column

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

# *** Read filtered VST matrix ***
vst_path <- "Perlo2022_counts_filters_VST_CNC_CV_above0.6.txt"
# Define first column as index
vst <- read.table(vst_path, header = TRUE, row.names = 1, check.names = FALSE) #encoding = "UTF-8", check.names = FALSE

### prolema está em como eu transformo o colnames do vst - evitar para ver como fica
# Extract only the first name before the hyphen in row names
#colnames(vst) <- sub("^([^_]+)_.*", "\\1", colnames(vst))
# This regular expression pattern looks for a lowercase letter followed by an uppercase letter and inserts a space between them
#colnames(vst) <- sub("([a-z])([A-Z])", "\\1 \\2", colnames(vst))
###

# *** Read metadata ***
metadata_path <-"infos_perlo_metadata.tsv"
metadata <- read.table(metadata_path, sep = "\t", header = T)

# *** Make vectors for each interesting module ***
for (i in Nmods$'Mod No'){
  assign(paste0("Module", i), modules %>% filter(module_No == i))
}

# *** Define groups to plot (Groups, Genotypes) ***
sample_table <- metadata

sample_table$Groups <- sub("^.*?_(\\S+)_\\d+-weeks.*$", "\\1", metadata$Run)
sample_table$Groups <- as.factor((sample_table$Groups))

sample_table$Replicate <- sub("^.*?_([^_]+_\\d+-weeks)_\\S+$", "\\1", metadata$Run)
sample_table$Replicate <- as.factor((sample_table$Replicate))

sample_table$Genotypes <- sub("^.*_(\\S+)$", "\\1", metadata$Run)
sample_table$Genotypes <- as.factor((sample_table$Genotypes))

# This regular expression pattern looks for a lowercase letter followed by an uppercase letter and inserts a space between them
#sample_table$Genotypes <- sub("([a-z])([A-Z])", "\\1 \\2", sample_table$Genotypes)

# Group (Genotypes and Groups)
#sample_table$Group <- as.factor(paste(sample_table$Genotypes, ' ', sample_table$Groups, ' ', sample_table$Replicate, sep=''))
#sample_table$Group <- as.factor(paste(sample_table$Groups, ' ', sample_table$Replicate, sep=''))
sample_table$Group <- as.factor(paste(sample_table$Replicate, sep=''))

annotation_col <- sample_table

# Remove duplicates in annotation_col (metadata)
unique_annotation_col <- distinct(annotation_col, Group, Genotypes, Groups, Replicate)
anot <- select(unique_annotation_col, Genotypes, Group)

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

# *** Reorder the rows of 'anot' based on the order of 'Genotypes' ***
merged_df <- merge(anot, data.frame(Genotypes = colnames(df)), by = "Genotypes", all.x = TRUE)

anot <- merged_df[order(match(colnames(df), merged_df$Genotypes)), ]

# Force the use of "-" in the rownames instead of "."
colnames(df) <- gsub("Internode_", "", colnames(df))
colnames(df) <- gsub("-weeks_", "_", colnames(df))
colnames(df) <- gsub("\\.", "-", colnames(df))

# Update 'anot' rownames -> Genotypes + Groups in names
rownames(anot) <- colnames(df)

# red (-) black (0) green (+)
my_palette = colorRampPalette(c("red", "black", "green"))(n=1000)

# *** Plot heatmap with mean values for column (i.e. mean for each module)
#png("meanOf40ModulesHeatmap.png", res = 300, width = 10*800, height = 10*2850) # Too big
png("meanOf40ModulesHeatmap.png", res = 300, width = 10*1500, height = 5*800)

pheatmap(df,
         main ="Genotypes contrasting in biomass production (CNC Modules)" ,
         scale = "row",
         annotation_col = anot,
         show_rownames = T,
         col = my_palette,
         cluster_cols = T,
         cluster_rows = T,
         cellwidth = NA,
         cellheight = 8,
         angle_col = 0)
dev.off()
