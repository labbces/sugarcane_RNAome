library(tidyverse)
library(DESeq2)
library(pheatmap)
library(viridis)
library(wesanderson)
library(RColorBrewer)

rm(list=ls())

DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/moduleHeatmap/Perlo/CNC"
setwd(DIR)

Nmods <- read.table("all_mods.txt", header = F)                                    # read file with modules numbers (1,2,3,4,5 ...) 
colnames(Nmods) <- "Mod No"

modules_path <- "Perlo2022_counts_filters_VST_topCV_mcl_formated_cliques.csv"      # read formated modules (cliques)

modules <- read.table(modules_path, header = F)                                    # remove duplicates from the first column
modules <- modules[!duplicated(modules$V1), ]

row.names(modules) <- modules$V1                                                   # set rownames as first column
modules <- modules[, -1]                                                           # remove first column
modules <- modules[, c(2,1)]                                                       # reorder columns

colnames(modules) <- c("module_No", "classification")
modules$gene <- rownames(modules)

vst_path <- "Perlo2022_counts_filters_VST_CNC_CV_above0.6.txt"                     # read filtered VST matrix
vst <- read.table(vst_path, header = TRUE, row.names = 1, check.names = FALSE)     # define first column as index

metadata_path <-"infos_perlo_metadata.tsv"                                         # read metadata 
metadata <- read.table(metadata_path, sep = "\t", header = T)

sample_table <- metadata                                                           # define groups to plot (Groups, Genotypes)

sample_table$Internode <- sub("^.*?_(\\S+)_\\d+-weeks.*$", "\\1", metadata$Run)    # find weeks in the columns
sample_table$Internode <- as.factor((sample_table$Internode))

sample_table$Replicate <- sub("^.*?_(\\d+)-weeks_\\S+$", "\\1", metadata$Run)      # week number (19 or 37)
sample_table$Replicate <- as.factor((sample_table$Replicate))

sample_table$Genotypes <- sub("^.*_(\\S+)$", "\\1", metadata$Run)                  # genotype (last name after underscore)
sample_table$Genotypes <- as.factor((sample_table$Genotypes))

sample_table$Group <- as.factor(paste(sample_table$Replicate, sep=''))             # group (Genotypes and Groups)

annotation_col <- sample_table

unique_annotation_col <- distinct(annotation_col, Genotypes, Internode, Replicate) # remove duplicates in annotation_col (metadata)

anot <- select(unique_annotation_col, Genotypes, Internode, Replicate)             # annotation columns

my_palette = colorRampPalette(c("red", "black", "green"))(n=1000)                  # red (-) black (0) green (+)

for (i in Nmods[,1]){
  names <- modules[modules$module_No == i,]
  df <- vst[names$gene,]
  
  # Reorder the rows of 'anot' based on the order of 'Genotypes'
  merged_df <- merge(anot, data.frame(Genotypes = colnames(df)), by = "Genotypes", all.x = TRUE)
  anot <- merged_df[order(match(colnames(df), merged_df$Genotypes)), ]
  
  # Rename colnames(df)
  colnames(df) <- gsub("Internode_([0-9]+)_([0-9]+).weeks_(.*)", "\\1_\\2_\\3", colnames(df))
  colnames(df) <- gsub("Internode_Ex\\.5_37\\.weeks_(.*)", "Ex.5_37_\\1", colnames(df))
  # Force the use of "-" in the rownames instead of "."
  colnames(df) <- gsub("\\.", "-", colnames(df))
  
  # update 'anot' rownames -> Genotypes + Groups in names
  rownames(anot) <- colnames(df)

  # pheatmap with mean values for column (mean condition expression) 
  png(paste0("module_", i, "_heatmap",".png", sep = ""), res = 300, width = 5*1500, height = 5*2850)
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