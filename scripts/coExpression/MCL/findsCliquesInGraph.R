#install.packages("igraph")

rm(list = ls())

library(igraph)

HOME_DIR = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/MCL"
setwd(HOME_DIR)

# cluster similarities
#data <- read.table("../fiberAndSugar/correr/MCL/CNC/clusterSimilarity/Correr_clusteringSimilarities.tsv", header = TRUE, sep = "\t")
#data <- read.table("../fiberAndSugar/hoang/MCL/CNC/clusterSimilarity/Hoang_clusteringSimilarities.tsv", header = TRUE, sep = "\t")
data <- read.table("../fiberAndSugar/perlo/MCL/CNC/clusterSimilarity/Perlo_clusteringSimilarities.tsv", header = TRUE, sep = "\t")

# *** Jaccard coefficient ***

# create Jaccard graph
graph_jaccard <- graph_from_data_frame(data, directed = FALSE, vertices = NULL)

# define edge weight as Jaccard coefficient
E(graph_jaccard)$weight <- data$Jaccard

# print graph
print(graph_jaccard)

# find cliques
jaccard_cliques <- weighted_cliques(graph_jaccard, vertex.weights = NULL, min.weight = 0.5, maximal = TRUE)
jaccard_cliques

# find largest cliques
jaccard_largest_cliques <- largest_weighted_cliques(graph_jaccard)
jaccard_largest_cliques

# print largest cliques size
weighted_clique_num(graph_jaccard)

# *** Overlap coefficient ***

# # create Overlap graph
graph_overlap <- graph_from_data_frame(data, directed = FALSE, vertices = NULL)

# define edge weight as Overlap coefficient
E(graph_overlap)$weight <- data$Overlap

# print graph
print(graph_overlap)

# find cliques
overlap_cliques <- weighted_cliques(graph_overlap, vertex.weights = NULL, min.weight = 0.5, maximal = TRUE)
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
save_cliques(jaccard_cliques, "cliques_jaccard.txt")
save_cliques(jaccard_largest_cliques, "largest_cliques_jaccard.txt")

# saving overlap cliques
save_cliques(overlap_cliques, "cliques_overlap.txt")
save_cliques(overlap_largest_cliques, "largest_cliques_overlap.txt")
