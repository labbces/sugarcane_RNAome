library(igraph, lib.loc="/Storage/data1/jorge.munoz/NRGSC.new/libraries")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript network_analysis_degree.r network_file")
}

# Read input
network_file <- args[1]

adjacency_list <- read.table(network_file, header = FALSE)
colnames(adjacency_list) <- c("V1", "V2", "weight")

G <- simplify(graph_from_data_frame(adjacency_list, directed = FALSE), remove.multiple = TRUE)

# Calculating and writing degree
names <- V(G)[order(degree(G), decreasing = TRUE)]
hub <- degree(G)[names]
names(hub) <- names(degree(G)[names])
hub <- as.data.frame(hub)

degree_output <- paste0(tools::file_path_sans_ext(network_file), "_degree.tsv")

write.table(rownames(hub), degree_output, col.names = FALSE, quote = FALSE, row.names = FALSE)
