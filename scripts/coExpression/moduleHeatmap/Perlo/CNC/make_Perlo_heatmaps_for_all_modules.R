library(tidyverse)
library(DESeq2)
library(pheatmap)
library(viridis)
library(wesanderson)
library(RColorBrewer)

rm(list=ls())

DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/moduleHeatmap/Perlo/CNC"
setwd(DIR)

#TODO: change to filtered_modules.txt
Nmods <- read.table("all_mods.txt", header = F)                                    # read file with modules numbers (1,2,3,4,5 ...) 
colnames(Nmods) <- "Mod No"

modules_path <- "Perlo2022_counts_filters_VST_topCV_mcl_formated_cliques.csv"      # read formated modules (cliques)

modules <- read.table(modules_path, header = F)                                    # read file with modules classification
colnames(modules) <- c("Gene", "Membership", "module_No")                          # rename columns

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

annotation_col <- sample_table

unique_annotation_col <- distinct(annotation_col, Genotypes, Internode, Replicate, Run) # remove duplicates in annotation_col (metadata)

anot <- select(unique_annotation_col, Genotypes, Internode, Replicate, Run)             # annotation columns

# make annotation vectors for each module 
for (i in Nmods$'Mod No'){
  assign(paste0("Module", i), modules %>% filter(module_No == i))
}

# make expression vectors for each module 
for (i in Nmods$'Mod No'){
  current_module <- modules %>% filter(module_No == i)
  rownames(current_module) <- current_module$Gene
  current_module <- current_module[, -1]      
  
  module_indices <- rownames(vst) %in% rownames(current_module)
  assign(paste0("module", i, "_dat"), vst[module_indices, ])
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

colnames(df) <- gsub("\\.", "-", colnames(df))                                     # force the use of "-" in the rownames instead of "."

# reorder the rows of 'anot' based on the order of 'Genotypes'
merged_df <- merge(data.frame(Run = colnames(df)), anot, by = "Run", all.x = TRUE)

anot <- merged_df[order(match(colnames(df), merged_df$Genotypes)), ]

#colnames(df) <- gsub("Internode_([0-9]+)_([0-9]+).weeks_(.*)", "\\1_\\2_\\3", colnames(df))  # rename colnames(df)
#colnames(df) <- gsub("Internode_Ex\\.5_37\\.weeks_(.*)", "Ex.5_37_\\1", colnames(df))

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