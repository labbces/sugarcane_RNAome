from random import sample
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

# Gold standard data of positive gene functional associations
# from https://www.inetbio.org/wormnet/downloadnetwork.php
#G = nx.read_edgelist("scripts/coExpression/networkMetrics/example/WormNet.v3.benchmark.txt")
#print(G)

# testing with my edgelist (cliques - similarity > 0.5)
#edge_list = "scripts/coExpression/networkMetrics/example/Correr_clusteringSimilarities_0.5.tsv"
edge_list = "scripts/coExpression/networkMetrics/example/Correr_clusteringSimilarities.tsv"

df = pd.read_csv(edge_list, sep='\t', usecols=[0,1,3])
print(df)
G = nx.from_pandas_edgelist(df, source='cluster_1', target='cluster_2', edge_attr='Overlap')
print(G)

# remove randomly selected nodes (to make example fast)
num_to_remove = int(len(G) / 1.5)
nodes = sample(list(G.nodes), num_to_remove)
G.remove_nodes_from(nodes)

# remove low-degree nodes
low_degree = [n for n, d in G.degree() if d < 10]
G.remove_nodes_from(low_degree)

# largest connected component
components = nx.connected_components(G)
largest_component = max(components, key=len)
H = G.subgraph(largest_component)

# compute centrality
centrality = nx.betweenness_centrality(H, k=10, endpoints=True) # k=10

# compute community structure
lpc = nx.community.label_propagation_communities(H)
community_index = {n: i for i, com in enumerate(lpc) for n in com}

#### draw graph ####
fig, ax = plt.subplots(figsize=(20, 15))
pos = nx.spring_layout(H, k=0.15, seed=4572321)
node_color = [community_index[n] for n in H]
node_size = [v * 20000 for v in centrality.values()]
nx.draw_networkx(
    H,
    pos=pos,
    with_labels=False,
    node_color=node_color,
    node_size=node_size,
    edge_color="gainsboro",
    alpha=0.4,
)

# Title/legend
font = {"color": "k", "fontweight": "bold", "fontsize": 10}
ax.set_title("Cliques association network (Correr2020 - Overlap coefficient)", font)
# Change font color for legend
font["color"] = "r"

ax.text(
    0.80,
    0.10,
    "node color = community structure",
    horizontalalignment="center",
    transform=ax.transAxes,
    fontdict=font,
)
ax.text(
    0.80,
    0.06,
    "node size = betweenness centrality",
    horizontalalignment="center",
    transform=ax.transAxes,
    fontdict=font,
)

# Resize figure for label readability
ax.margins(0.1, 0.05)
fig.tight_layout()
plt.axis("off")
plt.show()