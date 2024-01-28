library(igraph, lib.loc="/Storage/data1/jorge.munoz/NRGSC.new/libraries")

network_file="Perlo2022_mcl_out.txt"

adjacency_list <- read.table(network_file, header = F)
colnames(adjacency_list) <- c("V1", "V2", "weight")

G <- simplify(graph_from_data_frame(adjacency_list, directed = F), remove.multiple = T)

names <- V(G)[order(betweenness(G), decreasing = TRUE)]
hub <- betweenness(G)[names]
names(hub) <- names(betweenness(G)[names])
hub <- as.data.frame(hub)
write.table(rownames(hub), "Perlo2022_mcl_out_betweenness.tsv", col.names = F, quote = F, row.names = F)
