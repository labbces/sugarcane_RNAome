library(tidyverse)
library(DESeq2)
library(pheatmap)
library(viridis)
library(wesanderson)
library(RColorBrewer)

rm(list=ls())

DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/moduleHeatmap/Hoang/CNC"
setwd(DIR)

# *** Read file for modules numbers (1,2,3,4,5 ...) *** 
Nmods <- read.table("all_mods.txt", header = F)
#Nmods <- read.table("one_mods.txt", header = F)
colnames(Nmods) <- "Mod No"

# *** Read formated modules ***
# formated cliques
modules_path <- "Hoang2017_counts_filters_VST_topCV_mcl_formated_cliques.csv"

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
vst_path <- "Hoang2017_counts_filters_VST_CNC_CV_above1.2.txt"
# Define first column as index
vst <- read.table(vst_path, header = TRUE, row.names = 1, check.names = FALSE) #encoding = "UTF-8", check.names = FALSE

# Extract only the first name before the hyphen in row names
#colnames(vst) <- sub("^([^_]+)_.*", "\\1", colnames(vst))
# This regular expression pattern looks for a lowercase letter followed by an uppercase letter and inserts a space between them
#colnames(vst) <- sub("([a-z])([A-Z])", "\\1 \\2", colnames(vst))

# *** Read metadata ***
metadata_path <-"infos_hoang_metadata.tsv"
metadata <- read.table(metadata_path, sep = "\t", header = T, skip = 1)

# *** Define groups to plot (Groups, Genotypes) ***
sample_table <- metadata

# Groups
sample_table$Groups <- sub(".*_(top|bottom)-internode$", "\\1", metadata$Run)
sample_table$Groups <- as.factor((sample_table$Groups))

# Sugar content
sample_table$Sugar <- sub("^.*?_(.*?)_.*", "\\1", metadata$Run)
sample_table$Sugar <- as.factor((sample_table$Sugar))

# Genotypes
sample_table$Genotypes <- sub("^(.*?)_.*", "\\1", metadata$Run)
sample_table$Genotypes <- as.factor((sample_table$Genotypes)) 

# This regular expression pattern looks for a lowercase letter followed by an uppercase letter and inserts a space between them
#sample_table$Genotypes <- sub("([a-z])([A-Z])", "\\1 \\2", sample_table$Genotypes)

# Group (Genotypes and Groups)
sample_table$Group <- as.factor(paste(sample_table$Genotypes, ' ', sample_table$Groups, ' ', sample_table$Sugar, sep=''))

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




####


library(tidyverse)
library(DESeq2)
library(pheatmap)
library(viridis)
library(wesanderson)
library(RColorBrewer)

rm(list=ls())

DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/moduleHeatmap/Hoang/CNC"
setwd(DIR)

Nmods <- read.table("all_mods.txt", header = F)                                    # read file with modules numbers (1,2,3,4,5 ...) 
colnames(Nmods) <- "Mod No"

modules_path <- "Hoang2017_counts_filters_VST_topCV_mcl_formated_cliques.csv"      # read formated modules (cliques)

modules <- read.table(modules_path, header = F)                                    # remove duplicates from the first column
modules <- modules[!duplicated(modules$V1), ]

row.names(modules) <- modules$V1                                                   # set rownames as first column
modules <- modules[, -1]                                                           # remove first column
modules <- modules[, c(2,1)]                                                       # reorder columns

colnames(modules) <- c("module_No", "classification")
modules$gene <- rownames(modules)

vst_path <- "Hoang2017_counts_filters_VST_CNC_CV_above1.2.txt"                     # read filtered VST matrix
vst <- read.table(vst_path, header = TRUE, row.names = 1, check.names = FALSE)     # define first column as index

metadata_path <-"infos_hoang_metadata.tsv"                                         # read metadata 
metadata <- read.table(metadata_path, sep = "\t", header = T, skip = 1)

sample_table <- metadata                                                           # define groups to plot (Groups, Genotypes)

sample_table$Internode <- sub(".*_(top|bottom)-internode$", "\\1", metadata$Run)   # find internode types (top or bottom)
sample_table$Internode <- as.factor((sample_table$Internode))

sample_table$Brix <- sub("^.*?_(.*?)_.*", "\\1", metadata$Run)                     # find sugar content (high or low)
sample_table$Brix <- as.factor((sample_table$Brix))

sample_table$Genotypes <- sub("^(.*?)_.*", "\\1", metadata$Run)                    # find genotype name
sample_table$Genotypes <- as.factor((sample_table$Genotypes)) 

#sample_table$Group <- as.factor(paste(sample_table$Replicate, sep=''))             # group (Genotypes and Groups)

annotation_col <- sample_table

unique_annotation_col <- distinct(annotation_col, Genotypes, Internode, Brix)      # remove duplicates in annotation_col (metadata)

anot <- select(unique_annotation_col, Genotypes, Internode, Brix)                  # annotation columns

my_palette = colorRampPalette(c("red", "black", "green"))(n=1000)                  # red (-) black (0) green (+)

# TODO: acho que merged_df nao está funcionando - verificar em Perlo
for (i in Nmods[,1]){
  names <- modules[modules$module_No == i,]
  df <- vst[names$gene,]
  
  # reorder the rows of 'anot' based on the order of 'Genotypes'
  merged_df <- merge(anot, data.frame(Genotypes = colnames(df)), by = "Genotypes", all.x = TRUE)
  anot <- merged_df[order(match(colnames(df), merged_df$Genotypes)), ]
  
  # force the use of "-" in the rownames instead of "."
  colnames(df) <- gsub("\\.", "-", colnames(df))
  
  # update 'anot' rownames -> Genotypes + Groups in names
  rownames(anot) <- colnames(df)

  # pheatmap with mean values for column (mean condition expression) 
  png(paste0("module_", i, "_heatmap",".png", sep = ""), res = 300, width = 5*800, height = 5*2850)
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
