#install.packages("igraph")

rm(list = ls())

library(igraph)

HOME_DIR = "D://GitHub/sugarcane_RNAome/scripts/coExpression/MCL"
HOME_DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/MCL"
setwd(HOME_DIR)

# cluster similarities
data <- read.table("../fiberAndSugar/correr/MCL/CNC/clusterSimilarity/Correr_clusteringSimilarities.tsv", header = TRUE, sep = "\t")
data <- read.table("/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/hoang/MCL/CNC/CV1.0/clusterSimilarity/Hoang_clusteringSimilarities.tsv", header = TRUE, sep = "\t")
data <- read.table("../fiberAndSugar/perlo/MCL/CNC/clusterSimilarity/Perlo_clusteringSimilarities.tsv", header = TRUE, sep = "\t")

# *** Jaccard coefficient ***
# Filter Jaccard >= 0.7
jaccard_filtered <- data[which(data$Jaccard>=0.7),]

# create Jaccard graph
graph_jaccard <- graph_from_data_frame(jaccard_filtered, directed = FALSE, vertices = NULL)

# define edge weight as Jaccard coefficient
E(graph_jaccard)$weight <- jaccard_filtered$Jaccard

# print graph
print(graph_jaccard)

# find cliques
jaccard_cliques <- weighted_cliques(graph_jaccard, vertex.weights = NULL, min.weight = 0.7, maximal = TRUE)
jaccard_cliques

# find largest cliques
jaccard_largest_cliques <- largest_weighted_cliques(graph_jaccard)
jaccard_largest_cliques

# print largest cliques size
weighted_clique_num(graph_jaccard)

# *** Overlap coefficient ***
# Filter Overlap >= 0.7
overlap_filtered <- data[which(data$Overlap>=0.7),]

# # create Overlap graph
graph_overlap <- graph_from_data_frame(overlap_filtered, directed = FALSE, vertices = NULL)

# define edge weight as Overlap coefficient
E(graph_overlap)$weight <- overlap_filtered$Overlap

# print graph
print(graph_overlap)

# find cliques
overlap_cliques <- weighted_cliques(graph_overlap, vertex.weights = NULL, min.weight = 0.7, maximal = TRUE)
overlap_cliques

# find largest cliques
overlap_largest_cliques <- largest_weighted_cliques(graph_overlap)
overlap_largest_cliques

# print largest cliques size
weighted_clique_num(graph_overlap)

# *** Saving results ***

# function to save cliques to a text file
save_cliques <- function(cliques, filename) {
  # open the file for writing
  file <- file(filename, "w")
  
  # write each clique to the file
  for (clique in cliques) {
    cat(paste(clique$name, collapse = "\t"), "\n", file = file)  # use clique$name to get the names of vertices in the clique
  }
  
  # close the file
  close(file)
}

# saving jaccard cliques
save_cliques(jaccard_cliques, "cliques_jaccard_0.7.txt")
save_cliques(jaccard_largest_cliques, "largest_cliques_jaccard_0.7.txt")

# saving overlap cliques
save_cliques(overlap_cliques, "cliques_overlap_0.7.txt")
save_cliques(overlap_largest_cliques, "largest_cliques_overlap_0.7.txt")

plot(data)
