library(igraph, lib.loc="/Storage/data1/jorge.munoz/NRGSC.new/libraries")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript network_analysis_betweenness.r network_file")
}

# Read input
network_file <- args[1]

adjacency_list <- read.table(network_file, header = FALSE)
colnames(adjacency_list) <- c("V1", "V2", "weight")

G <- simplify(graph_from_data_frame(adjacency_list, directed = FALSE), remove.multiple = TRUE)

# Calculating and writing betweenness
names <- V(G)[order(betweenness(G), decreasing = TRUE)]
hub <- betweenness(G)[names]
names(hub) <- names(betweenness(G)[names])
hub <- as.data.frame(hub)

betweenness_output <- paste0(tools::file_path_sans_ext(network_file), "_betweenness.tsv")

write.table(rownames(hub), betweenness_output, col.names = FALSE, quote = FALSE, row.names = FALSE)
