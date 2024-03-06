import networkx as nx
import plotly.graph_objects as go

file_path = '/home/felipe/Documents/ecoli-co-expression-network/results/cliques_overlap.txt'

G = nx.Graph()

with open(file_path, 'r') as file:
    for line in file:
        nodes = line.strip().split('\t')

        # If there is only one gene in the cluster, add it as an isolate node
        #if len(nodes) == 1:
            #G.add_node(nodes[0])
        #elif len(nodes) >= 2:    

        # If there is only one gene in the cluster, skip it
        if len(nodes) < 2:
            continue
        G.add_edges_from([(nodes[i], nodes[i + 1]) for i in range(len(nodes) - 1)])

# Get node positions with spring layout
pos = nx.spring_layout(G)

# Create interactive graph
fig = go.Figure()

# Add nodes and edges
for edge in G.edges():
    x0, y0 = pos[edge[0]]
    x1, y1 = pos[edge[1]]
    
    # Add edges
    fig.add_trace(go.Scatter(
        x=[x0, x1, None],
        y=[y0, y1, None],
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines'))

# Add isolate nodes
for node in G.nodes():
    x, y = pos[node]
    
    # Add isolate nodes
    fig.add_trace(go.Scatter(
        x=[x],
        y=[y],
        mode='markers',
        hoverinfo='text',
        marker=dict(size=10),
        text=node
    ))

# Add legend for all nodes
fig.add_trace(go.Scatter(
    x=[],
    y=[],
    mode='markers',
    hoverinfo='text',
    marker=dict(
        size=10,
        colorbar=dict(
            thickness=15,
            title='Node Connections',
            xanchor='left',
            titleside='right'
        )
    ),
    showlegend=True 
))

# Update layout
fig.update_layout(
    hovermode='closest',
    margin=dict(b=0, l=0, r=0, t=0),
    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))

fig.show()