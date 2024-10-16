# read modules
Nmods <- read.table("Correr2020_counts_filters_VST_topCV_mcl_formated_cliques_filtered.csv", header = F)
colnames(Nmods) <- "Mod No"

# files
modules_path <- "Correr2020_counts_filters_VST_topCV_mcl_formated_cliques.csv"
vst_path <- "Correr2020_counts_filters_VST_CNC_CV_above_2.txt"
metadata_path <-"metadata.tsv"

# libraries
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(viridis)
library(wesanderson)

#modules <- read.table(modules_path, row.names = 1, header = F)
# remove duplicated genes
modules <- read.table(modules_path, header = FALSE)
modules <- modules[!duplicated(modules[, 1]), ]
row.names(modules) <- modules[, 1]
modules <- modules[, -1]

colnames(modules) <- c("pertinence", "module_No")

# read full vst matrix 
vst <- read.table(vst_path, row.names = 1, header = T, check.names = FALSE)
metadata <- read.table(metadata_path, sep = "\t", header = T)

# Make vectors for each intersting module
for (i in Nmods$'Mod No'){
assign(paste0("Module", i), modules %>% filter(module_No == i))
}

sample_table <- metadata
sample_table$Trait <- as.factor((sample_table$Trait))
sample_table$Genotype <- as.factor((sample_table$Genotype))
sample_table$Group <- as.factor(paste(sample_table$Genotype, '_', sample_table$Trait, sep=''))

#df <- as.data.frame((assay(vst)))
annotation_col <- sample_table
anot <- select(annotation_col, "Group", "Genotype", "Trait")

for (i in Nmods$'Mod No'){
assign(paste0("module", i, "_dat"), vst[rownames(eval(as.name(paste0("Module",i)))),])
}

trait_genotype <- as.numeric(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
trait_biomass_condition <- as.numeric(c("0", "1", "1", "1", "1", "0", "0", "0", "0", "0", "1", "0"))
trait_interaction <- trait_genotype * trait_biomass_condition

trait_fiber <- as.numeric(c("9.83", "25.55", "24.64", "22.35", "22.22", "11.14", "12", "21.86", "12.14", "12.97", "18.96", "10.82"))
trait_brix <- as.numeric(c("16.53", "14.54", "14.78", "15.91", "12.03", "20.52", "21.01", "13.99", "21.29", "22.76", "17.61", "19.6"))

# get_cor_module(module10_dat, trait_genotype, "M10")
get_cor_module <- function(x,trait, module){
df <- x
#df2 <- data.frame(matrix(ncol = 12, nrow = dim(df)[1]))
#colnames(df2) <- unique(anot$Group)
#rownames(df2) <- rownames(df)
# removed lines to colapse replicates in original code
pca <- prcomp(t(df))
eigen <-pca$x[,1]
cor <- cor.test(eigen, trait, method = "spearman", exact=FALSE)
pvalue <- cor$p.value
rho <- cor$estimate

df_out <- data.frame(rho=rho,
p = pvalue,
Module = module,
row.names = NULL
)
return(df_out)
}

condition_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))
genotype_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))

fiber_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))
brix_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))
interaction_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))

colnames(condition_table) <- c("rho", "pvalue", "module")
colnames(genotype_table) <- c("rho", "pvalue", "module")

colnames(fiber_table) <- c("rho", "pvalue", "module")
colnames(brix_table) <- c("rho", "pvalue", "module")
colnames(interaction_table) <- c("rho", "pvalue", "module")

for (i in Nmods$'Mod No'){
condition_table[i,] <- get_cor_module(eval(as.name(paste0("module",i,"_dat"))), trait_biomass_condition, paste0("M", i))

genotype_table[i,] <- get_cor_module(eval(as.name(paste0("module",i,"_dat"))), trait_genotype, paste0("M", i))

fiber_table[i,] <- get_cor_module(eval(as.name(paste0("module",i,"_dat"))), trait_fiber, paste0("M", i))
brix_table[i,] <- get_cor_module(eval(as.name(paste0("module",i,"_dat"))), trait_brix, paste0("M", i))
interaction_table[i,] <- get_cor_module(eval(as.name(paste0("module",i,"_dat"))), trait_interaction, paste0("M", i))
}
 
filtered_genotype <- mutate(genotype_table, adjusted = p.adjust(genotype_table$pvalue, method = "bonferroni")) %>% mutate(genotype_table, fdr = p.adjust(genotype_table$pvalue, method = "fdr")) %>% filter(adjusted  < 0.05 & abs(rho) > 0.07)

filtered_condition <- mutate(condition_table, adjusted = p.adjust(condition_table$pvalue, method = "bonferroni")) %>% mutate(condition_table, fdr = p.adjust(condition_table$pvalue, method = "fdr")) %>% filter(adjusted  < 0.05 & abs(rho) > 0.7)

filtered_interaction <- mutate(interaction_table, adjusted = p.adjust(interaction_table$pvalue, method = "bonferroni")) %>% mutate(interaction_table, fdr = p.adjust(interaction_table$pvalue, method = "fdr")) %>% filter(adjusted < 0.05 & abs(rho) > 0.7)

filtered_fiber <- mutate(fiber_table, adjusted = p.adjust(fiber_table$pvalue, method = "bonferroni")) %>% mutate(fiber_table, fdr = p.adjust(fiber_table$pvalue, method = "fdr")) %>% filter(adjusted < 0.05 & abs(rho) > 0.7)

filtered_brix <- mutate(brix_table, adjusted = p.adjust(brix_table$pvalue, method = "bonferroni")) %>% mutate(brix_table, fdr = p.adjust(brix_table$pvalue, method = "fdr")) %>% filter(adjusted < 0.05 & abs(rho) > 0.7)

write.table(filtered_genotype, "genotype_bonferroni_005_rho_07.txt", row.names = F, quote = F, col.names = T)
write.table(filtered_condition, "condition_bonferroni_005_rho_07.txt", row.names = F, quote = F, col.names = T)
write.table(filtered_interaction, "interaction_bonferroni_005_rho_07.txt", row.names = F, quote = F, col.names = T)
write.table(filtered_fiber, "fiber_bonferroni_005_rho_07.txt", row.names = F, quote = F, col.names = T)
write.table(filtered_brix, "brix_bonferroni_005_rho_07.txt", row.names = F, quote = F, col.names = T)

library(ggVennDiagram)
library(ggplot2)

df1 <- read.table("genotype_bonferroni_005_rho_07.txt", header = True)
df2 <- read.table("condition_bonferroni_005_rho_07.txt", header = True)
df3 <- read.table("interaction_bonferroni_005_rho_07.txt", header = True)
df4 <- read.table("fiber_bonferroni_005_rho_07.txt", header = True)
df5 <- read.table("brix_bonferroni_005_rho_07.txt", header = True)

all <- list(genotype = df1$module, replicate = df2$module, internode = df3$Module, interaction = df4$module, df5$module)
all_upset_plot <- ggVennDiagram(all, force_upset = T ,order.set.by = "none", nintersects = 17)
ggsave("upsetplot_correlations.png", all_upset_plot, "png", dpi = 300)
