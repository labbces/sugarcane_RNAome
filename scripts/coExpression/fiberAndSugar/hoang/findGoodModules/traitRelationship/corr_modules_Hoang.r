# read modules
Nmods <- read.table("Hoang2017_counts_filters_VST_topCV_mcl_formated_cliques_filtered.csv", header = F)
colnames(Nmods) <- "Mod No"

# files
modules_path <- "Hoang2017_counts_filters_VST_topCV_mcl_formated_cliques.csv"
vst_path <- "Hoang2017_counts_filters_VST_CNC_CV_above1.2.txt"
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
sample_table$Internode <- as.factor((sample_table$Internode))
sample_table$Genotype <- as.factor((sample_table$Genotype))
sample_table$Condition <- as.factor((sample_table$Condition))
sample_table$Group <- as.factor(paste(sample_table$Genotype, '_', sample_table$Condition, "_", sample_table$Internode, sep=''))

#df <- as.data.frame((assay(vst)))
annotation_col <- sample_table
anot <- select(annotation_col, "Group", "Genotype", "Internode", "Condition")

for (i in Nmods$'Mod No'){
assign(paste0("module", i, "_dat"), vst[rownames(eval(as.name(paste0("Module",i)))),])
}

# transform metadata to numeric
data <- metadata
data_numeric <- data %>%
mutate(across(everything(), ~ as.integer(as.factor(.))))

trait_genotype <- as.numeric(data_numeric$Genotype)
trait_condition <- as.numeric(data_numeric$Condition)
trait_internode <- as.numeric(data_numeric$Internode)
trait_interaction <- trait_genotype * trait_internode

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

genotype_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))
condition_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))
internode_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))
interaction_table <- data.frame(matrix(ncol = 3, nrow = dim(Nmods)[1]))

colnames(genotype_table) <- c("rho", "pvalue", "module")
colnames(condition_table) <- c("rho", "pvalue", "module")
colnames(internode_table) <- c("rho", "pvalue", "module")
colnames(interaction_table) <- c("rho", "pvalue", "module")

for (i in Nmods$'Mod No') {
module_name <- paste0("module", i, "_dat")
module_data <- eval(as.name(module_name))
  
if (nrow(module_data) > 5) {
genotype_table[i,] <- get_cor_module(module_data, trait_genotype, paste0("M", i))
condition_table[i,] <- get_cor_module(module_data, trait_condition, paste0("M", i))
internode_table[i,] <- get_cor_module(module_data, trait_internode, paste0("M", i))
interaction_table[i,] <- get_cor_module(module_data, trait_interaction, paste0("M", i))
} else {
message(paste("Module", i, "has less than 5 genes and was ignored."))
}
}

filtered_genotype <- mutate(genotype_table, adjusted = p.adjust(genotype_table$pvalue, method = "bonferroni")) %>% mutate(genotype_table, fdr = p.adjust(genotype_table$pvalue, method = "fdr")) %>% filter(adjusted  < 0.05 & abs(rho) > 0.1)

filtered_condition <- mutate(condition_table, adjusted = p.adjust(condition_table$pvalue, method = "bonferroni")) %>% mutate(condition_table, fdr = p.adjust(condition_table$pvalue, method = "fdr")) %>% filter(adjusted  < 0.05 & abs(rho) > 0.1)

filtered_internode <- mutate(internode_table, adjusted = p.adjust(internode_table$pvalue, method = "bonferroni")) %>% mutate(internode_table, fdr = p.adjust(internode_table$pvalue, method = "fdr")) %>% filter(adjusted  < 0.05 & abs(rho) > 0.1)

filtered_interaction <- mutate(interaction_table, adjusted = p.adjust(interaction_table$pvalue, method = "bonferroni")) %>% mutate(interaction_table, fdr = p.adjust(interaction_table$pvalue, method = "fdr")) %>% filter(adjusted  < 0.05 & abs(rho) > 0.1)

#filtered <- filter(table, pvalue < 0.01 & abs(rho) > 0.7)
#write.table(filtered_condition, "../../results/perNCONDITION/module_nitrogen_condition_correlation_fdr_0.5.csv", row.names = F, quote = F) 
#write.table(filtered_genotype, "../../results/perNCONDITION/module_genotype_correlation_fdr_0.5.csv", row.names = F, quote = F)
#write.table(filtered_leaf_segment, "../../results/perNCONDITION/module_leaf_segment_fdr_0.5.csv", row.names = F, quote = F)
# Write tables

write.table(filtered_genotype, "genotype_bonferroni_005_rho_06.txt", row.names = F, quote = F, col.names = T)
write.table(filtered_week, "week_bonferroni_005_rho_06.txt", row.names = F, quote = F, col.names = T)
write.table(filtered_internode, "internode_bonferroni_005_rho_06.txt", row.names = F, quote = F, col.names = T)
write.table(filtered_interaction, "interaction_bonferroni_005_rho_06.txt", row.names = F, quote = F, col.names = T)

library(ggVennDiagram)
library(ggplot2)

df1 <- read.table("genotype_bonferroni_005_rho_06.txt", header = True)
df2 <- read.table("week_bonferroni_005_rho_06.txt", header = True)
df3 <- read.table("internode_bonferroni_005_rho_06.txt", header = True)
df4 <- read.table("interaction_bonferroni_005_rho_06.txt", header = True)

all <- list(genotype = df1$module, replicate = df2$module, internode = df3$Module, interaction = df4$module)
all_upset_plot <- ggVennDiagram(all, force_upset = T ,order.set.by = "none", nintersects = 17)
ggsave("upsetplot_correlations.png", all_upset_plot, "png", dpi = 300)
