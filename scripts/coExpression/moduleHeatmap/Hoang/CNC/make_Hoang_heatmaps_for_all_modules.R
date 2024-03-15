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

modules_path <- "Hoang2017_counts_filters_VST_topCV_mcl_formated_cliques.csv"      # read formated modules (cliques)

modules <- read.table(modules_path, header = F)                                    # remove duplicates from the first column
modules <- modules[!duplicated(modules$V1), ]

row.names(modules) <- modules$V1                                                   # set rownames as first column
modules <- modules[, -1]                                                           # remove first column
modules <- modules[, c(2,1)]                                                       # reorder columns

colnames(modules) <- c("module_No", "classification")

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

annotation_col <- sample_table

unique_annotation_col <- distinct(annotation_col, Genotypes, Internode, Brix)      # remove duplicates in annotation_col (metadata)

anot <- select(unique_annotation_col, Genotypes, Internode, Brix)                  # annotation columns

# make annotation vectors for each module 
for (i in Nmods$'Mod No'){
  assign(paste0("Module", i), modules %>% filter(module_No == i))
}

# make expression vectors for each module 
for (i in Nmods$'Mod No'){
  assign(paste0("module", i, "_dat"), vst[rownames(eval(as.name(paste0("Module",i)))),])
}

# calculate mean expression for each module
for (i in Nmods$'Mod No'){
  assign(paste0("dat", i),as.data.frame(t(colMeans(eval(as.name(paste0("module",i, "_dat"))))))) 
  labelling <-get(paste0("dat", i, sep = ""))
  labels <- paste0("Module ", i, sep = "")
  rownames(labelling) <- labels
  assign(paste0("dat", i, sep= ""), labelling)
}

dat_list <- list()                                                                 # create a list to aggregate all modules 

# make mean expression vectors for each module 
for (i in Nmods$'Mod No'){
  dat_list[[i]] <- eval(as.name(paste0("dat",i)))
}

df <- do.call("rbind", dat_list)                                                   # create a DataFrame from all modules list

# reorder the rows of 'anot' based on the order of 'Genotypes'
merged_df <- merge(anot, data.frame(Genotypes = colnames(df)), by = "Genotypes", all.x = TRUE) 

anot <- merged_df[order(match(colnames(df), merged_df$Genotypes)), ]

colnames(df) <- gsub("\\.", "-", colnames(df))                                     # force the use of "-" in the rownames instead of "."

rownames(anot) <- colnames(df)                                                     # update 'anot' rownames -> Genotypes + Groups in names

my_palette = colorRampPalette(c("red", "black", "green"))(n=1000)                  # red (-) black (0) green (+)

# pheatmap with mean values for each module (mean module expression) 
png("meanOf3ModulesHeatmap.png", res = 300, width = 10*1500, height = 5*800)

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