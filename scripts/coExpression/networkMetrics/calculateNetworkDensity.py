'''
A densidade de uma rede (ou grafo) é uma medida que indica o quão próximo a rede está de ser um grafo completo (onde todos os nós estão conectados entre si). 
A densidade D de uma rede não direcionada pode ser calculada usando a fórmula:
  
Densidade = 2E / N(N - 1)

 Onde: 
- E é o número de arestas.
- N é o número de nós.
'''

def calculate_density(nodes, edges):
    return (2 * edges) / (nodes * (nodes - 1))

networks = {
    "Hoang_CV_1.0": {"nodes": 594573, "edges": 6972653560},
    "Hoang_CV_1.2": {"nodes": 240318, "edges": 1698491712},
    "Correr_CV_1.5": {"nodes": 992460, "edges": 33696299438},
    "Correr_CV_2.0": {"nodes": 255808, "edges": 3443744874},
    "Perlo_CV_0.6": {"nodes": 304865, "edges": 76956069},
}

densities = {name: calculate_density(data["nodes"], data["edges"]) for name, data in networks.items()}
print(densities)