library(tidyverse)
library(DESeq2)
library(pheatmap)
library(viridis)
library(wesanderson)
library(RColorBrewer)

rm(list=ls())

DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/moduleHeatmap/Hoang/CNC"
setwd(DIR)

#TODO: change to filtered_modules.txt
Nmods <- read.table("all_mods.txt", header = F)                                    # read file with modules numbers (1,2,3,4,5 ...) 
colnames(Nmods) <- "Mod No"

modules_path <- "Hoang2017_counts_filters_VST_topCV_mcl_formated_cliques.csv"      # filename for formated modules (cliques)
classification_path <- "Hoang2017_counts_filters_VST_CV_classified.txt"            # filename for annotation file (CV and Function)

modules <- read.table(modules_path, header = F)                                    # read file with modules classification
colnames(modules) <- c("Gene", "Membership", "module_No")                          # rename columns

classification_file <- read.table(classification_path, row.names = 1, header = T)  # read file with function and pantranscriptome classfication

merged_annotation <- merge(modules, classification_file, by.x = "Gene", by.y = "row.names", all.x = TRUE)

vst_path <- "Hoang2017_counts_filters_VST_CNC_CV_above1.2.txt"                     # read filtered VST matrix
vst <- read.table(vst_path, header = TRUE, row.names = 1, check.names = FALSE)     # define first column as index

metadata_path <-"infos_hoang_metadata.tsv"                                         # read metadata 
metadata <- read.table(metadata_path, sep = "\t", header = T, skip = 1)

sample_table <- metadata                                                           # define groups to plot (Groups, Genotypes)

sample_table$Internode <- sub(".*_(top|bottom)-internode$", "\\1", metadata$Run)   # find internode types (top or bottom)
sample_table$Internode <- as.factor((sample_table$Internode))

sample_table$Brix <- sub("^.*?_(.*?)_.*", "\\1", metadata$Run)                     # find sugar content (high or low)
sample_table$Brix <- gsub("\\-", " ", sample_table$Brix)
sample_table$Brix <- as.factor((sample_table$Brix))

sample_table$Genotypes <- sub("^(.*?)_.*", "\\1", metadata$Run)                    # find genotype name
sample_table$Genotypes <- as.factor((sample_table$Genotypes)) 

annotation_col <- sample_table

unique_annotation_col <- distinct(annotation_col, Genotypes, Internode, Brix, Run)      # remove duplicates in annotation_col (metadata)

anot <- select(unique_annotation_col, Genotypes, Internode, Brix, Run)                  # annotation columns

my_palette = colorRampPalette(c("red", "black", "green"))(n=1000)                       # red (-) black (0) green (+)

for (i in Nmods[,1]){
  names <- merged_annotation[merged_annotation$module_No == i,]
  df <- vst[names$Gene,]
  
  # force the use of "-" in the rownames instead of "."
  colnames(df) <- gsub("\\.", "-", colnames(df))
  
  # reorder the rows of 'anot' based on the order of 'Genotypes'
  merged_df <- merge(data.frame(Run = colnames(df)), anot, by = "Run", all.x = TRUE)
  anot_col <- merged_df[order(match(colnames(df), merged_df$Genotypes)), ]
  
  # create annotation row object
  anot_row <- names[, c("Gene", "Membership", "Classification", "Function")]
  rownames(anot_row) <- names$Gene
  anot_row <- anot_row[, -1]                                                       # remove first column
  colnames(anot_row) <- c("Membership", "Group", "Function")                       # rename row annotation columns  
  
  # update 'anot' rownames -> Genotypes + Groups in names
  rownames(anot_col) <- colnames(df)
  anot_col <- anot_col[, -1]                                                       # remove first column
  colnames(anot_col) <- c("Genotypes", "Internode", "Sugar Content")
  
  # pheatmap with mean values for column (mean condition expression) 
  png(paste0("module_", i, "_heatmap",".png", sep = ""), res = 300, width = 5*800, height = 5*2850)
  pheatmap(df,
           main =paste0("Genotypes contrasting in biomass production (CNC Module ",i, ")", sep = "") ,
           scale = "row",
           annotation_col = anot_col,
           annotation_row = anot_row,
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