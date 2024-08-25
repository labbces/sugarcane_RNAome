install.packages("ggVennDiagram")

library(ggVennDiagram)
library(ggplot2)

dir="C:/Users/PC/Desktop/temp/traitRelationship"
setwd(dir)

list.files()

df1 <- read.table("genotype_bonferroni_005_rho_06.txt", header = T)
df2 <- read.table("week_bonferroni_005_rho_06.txt", header = T)
df3 <- read.table("internode_bonferroni_005_rho_06.txt", header = T)
df4 <- read.table("interaction_bonferroni_005_rho_06.txt", header = T)

all <- list("genótipo" = df1$module, "estágio de desenvolvimento" = df2$module, "entrenó" = df3$module)#, interaction = df4$module)
all_upset_plot <- ggVennDiagram(all, force_upset = T ,order.set.by = "none", nintersects = 3,  top.bar.y.label = NULL,  sets.bar.x.label = "Módulos", intersection.matrix.color = "grey30")


#ggsave("upsetplot_correlations.png", all_upset_plot, "png", dpi = 300)

