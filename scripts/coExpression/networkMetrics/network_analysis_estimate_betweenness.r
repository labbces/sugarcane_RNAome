library(igraph, lib.loc="/Storage/data1/jorge.munoz/NRGSC.new/libraries")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript network_analysis_estimate_betweenness.r network_file")
}

# Read input
network_file <- args[1]

adjacency_list <- read.table(network_file, header = FALSE)
colnames(adjacency_list) <- c("V1", "V2", "weight")

G <- simplify(graph_from_data_frame(adjacency_list, directed = FALSE), remove.multiple = TRUE)

# Calculating and writing estimate betweenness
names <- V(G)[order(estimate_betweenness(G, cutoff = 2), decreasing = TRUE)]
hub <- estimate_betweenness(G, cutoff = 2)[names]
names(hub) <- names(estimate_betweenness(G, cutoff = 2)[names])
hub <- as.data.frame(hub)

betweenness_output <- paste0(tools::file_path_sans_ext(network_file), "_estimate_betweenness.tsv")

write.table(rownames(hub), betweenness_output, col.names = FALSE, quote = FALSE, row.names = FALSE)
