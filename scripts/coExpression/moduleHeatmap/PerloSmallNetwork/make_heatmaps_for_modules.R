library(tidyverse)
library(DESeq2)
library(pheatmap)
library(viridis)
library(wesanderson)

rm(list=ls())

DIR = "/home/felipe/Documents/exploring_Jorge_heatmaps/Perlo/CNC/heatmap_per_module"
setwd(DIR)

# *** files ***

# Read nitrogen modules 
# Este arquivo tem a quantidade de modulos (qtd de plots a serem criados)
# wc -l /Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/pearson_mcl/out.1.8.number 
# 199 (199 modules)

# Este primeiro teste tem apenas 4 (1,2,3,4)
Nmods <- read.table("all_mods.txt", header = F)
colnames(Nmods) <- "Mod No"

# gene module.number
modules_path <- "Perlo2022_mcl_out_I2.0.formated.csv"

# vst counts
# Nao estou usando VST neste teste, apenas raw TPM (grep genes do modules_path na matriz do Renato para Perlo)
vst_path <- "OnlyPerloGenes_MCLmodules_Perlo2022_tpm_matrix.txt"

# metadata (Sample Number,Sample Name,Sample_file,Condition,Genotype,DevStage,Individual,Group)
metadata_path <-"infos_perlo_metadata.tsv"

# *** read files ***

# read modules
modules <- read.table(file.path(DIR, modules_path), row.names = 1, header = F)
modules$gene <- rownames(modules)
colnames(modules) <- c("module_No")
  
# read full vst matrix 
vst <- read.table(vst_path)

# read metadata
metadata <- read.table(metadata_path, sep = "\t", header = T)
sample_table <- metadata

#sample_table$Condition <- as.factor((sample_table$Condition))
sample_table$Condition <- sub("^.*?_(\\S+)_\\d+-weeks.*$", "\\1", metadata$Run)
sample_table$Condition <- as.factor((sample_table$Condition))

#sample_table$Genotype <- as.factor((sample_table$Genotype))
sample_table$Genotype <- sub("^.*_(\\S+)$", "\\1", metadata$Run)
sample_table$Genotype <- as.factor((sample_table$Genotype)) 

#sample_table$DevStage <- as.factor((sample_table$DevStage))
sample_table$DevStage <- sub("^.*?_([^_]+_\\d+-weeks)_\\S+$", "\\1", metadata$Run)
sample_table$DevStage <- as.factor((sample_table$DevStage))

#sample_table$Individual <- as.factor((sample_table$Individual))

sample_table$Group <- as.factor(paste(sample_table$Genotype, '_', sample_table$Condition, '_', sample_table$DevStage, sep=''))

#df <- as.data.frame((assay(vst)))
annotation_col <- sample_table
anot <- select(annotation_col, "Group","DevStage","Genotype", "Condition")

# *** PLOT HEATMAP 

# Subset data based on module number
module_genes <- modules[modules$module_No == 1, "gene"]
df <- vst[vst[1]$V1 %in% module_genes, ]

# Set column names
colnames(df) <- colnames(vst)
  
# Set row names for annotation
#rownames(anot) <- colnames(df)
  
# Create a heatmap
png(paste0(DIR, "/hp_samples_module_", 1, ".png"), res = 300, width = 5.5 * 800, height = 1.5 * 2850)
pheatmap(df,
         main = paste0("Nitrogen Module ", i)
         #scale = "row",
         #annotation_col = anot,
         #show_rownames = TRUE,
         #col = inferno(299),
         #cluster_cols = TRUE,
         #cluster_rows = TRUE,
         #cellwidth = NA,
         #cellheight = 8) 
dev.off()

####



# codigo jorge
for (i in Nmods[,1]){
  df <- vst[modules[modules$module_No == i,]]
  names <- modules[modules$module_No == i,]
  df <- vst[names$gene,]
  colnames(df) <- colnames(vst)
  rownames(anot) <- colnames(df)
  # heat map with ean values for column i.e media por modulo
  png(paste0(DIR, "hp_samples_module_", i, ".png", sep = ""), res = 300, width = 5.5*800, height = 1.5*2850)
  pheatmap(df,
           main =paste0("Nitrogen Module ",i, sep = "") ,
           scale = "row",
           annotation_col = anot,
           show_rownames = T,
           col = inferno(299),
           cluster_cols = T,
           cluster_rows = T,
           cellwidth = NA,
           cellheight = 8) 
  dev.off()
  print(i)
}
  
  
  df2 <- data.frame(matrix(ncol = 16, nrow = dim(df)[1]))
  colnames(df2) <- unique(anot$Group)
  rownames(df2) <- rownames(df)
  
  df2$NR_10_B <- (df$Sample_38 + df$Sample_42 + df$Sample_46)/3
  df2$NR_10_B0 <- (df$Sample_37 + df$Sample_41 + df$Sample_45)/3
  df2$NR_10_M <- (df$Sample_39 + df$Sample_43 + df$Sample_47)/3
  df2$NR_10_P <- (df$Sample_40 + df$Sample_44 + df$Sample_48)/3
  df2$NR_270_B <- (df$Sample_26 + df$Sample_30 + df$Sample_34)/3
  df2$NR_270_B0 <- (df$Sample_25 + df$Sample_29 + df$Sample_33)/3
  df2$NR_270_M <- (df$Sample_27 + df$Sample_31 + df$Sample_35)/3
  df2$NR_270_P <- (df$Sample_28 + df$Sample_32 + df$Sample_36)/3
  df2$R_10_B <- (df$Sample_14 + df$Sample_18 + df$Sample_22)/3
  df2$R_10_B0 <- (df$Sample_13 + df$Sample_17 + df$Sample_21)/3
  df2$R_10_M  <- (df$Sample_15 + df$Sample_19 + df$Sample_23)/3
  df2$R_10_P <- (df$Sample_16 + df$Sample_20 + df$Sample_24)/3
  df2$R_270_B <- (df$Sample_2 + df$Sample_6 + df$Sample_10)/3
  df2$R_270_B0 <- (df$Sample_1 + df$Sample_5 + df$Sample_9)/3
  df2$R_270_M <- (df$Sample_3 + df$Sample_7 + df$Sample_11)/3
  df2$R_270_P <- (df$Sample_4 + df$Sample_8 + df$Sample_12)/3
  
  #write.table(df2, "/Storage/data1/jorge.munoz/NRGSC.new/networks/results/perNCONDITION/all_matrix_for_pca.csv", sep = ",")
  
  # metadata (Genotype, Condition, DevStage)
  anot2 <- read.table("/Storage/data1/jorge.munoz/NRGSC.new/networks/data/metadata_complete_summarized.csv", sep = ",", header = T, colClasses = c("factor", "factor", "factor" ), row.names = 1)
  
  colorDev=wes_palette("Zissou1", 4, type = "discrete")
  colorGen=wes_palette("Darjeeling1", 2, type = "discrete")
  colorCond=wes_palette("Chevalier1", 2, type = "discrete")
  
  ann_colors = list(
    DevStage = c(B=colorDev[1], B0=colorDev[2], M=colorDev[3], P=colorDev[4]),
    Genotype = c(R = colorGen[1], NR = colorGen[2]),
    Condition = c("10" = colorCond[1], "270" = colorCond[2]))
  
  png(paste0(DIR, "hp_grouped_module_" , i, ".png", sep = ""), res = 300, width = 2*800, height = 1.5*2850)
  pheatmap(df2, 
           main ="Nitrogen Modules",
           cellheight = 8,
           cellwidth =NA,
           scale = "row",
           show_rownames = T,
           col = inferno(299),
           cluster_cols = T,
           cluster_rows = T,
           annotation_col  = anot2,
           annotation_colors = ann_colors) 
  dev.off()
  print(i)
}
