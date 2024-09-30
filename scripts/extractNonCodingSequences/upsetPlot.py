#!/usr/bin/env python

import pandas as pd
import upsetplot
import matplotlib.pyplot as plt
import itertools
import numpy as np

with open('CPC2_ncRNAs_list.txt') as f:
    set1 = set(f.read().strip().split(','))

with open('CPC2_PLncPRO_ncRNAs_list.txt') as f:
    set2 = set(f.read().strip().split(','))

with open('RNAplonc_ncRNAs_list.txt') as f:
    set3 = set(f.read().strip().split(','))

# dicionario com as categorias
genes_dict = {"CPC2": set1, "PLncPRO": set2, "RNAplonc": set3}
categories = list(genes_dict.keys())

comb_list_list = []
comb_intersection_length_list = []

# display intersection: https://stackoverflow.com/questions/76230150/how-to-display-intersection-values-instead-of-distinct-values-in-upset-plot

# identificando o tamanho da intersecao para cada combinacao de categorias
for i in range(len(categories)):
    comb_list = list(itertools.combinations(categories, i + 1))
    for elem in comb_list:
        comb_list_list.append(elem)
        cat_lists = [genes_dict[x] for x in elem]
        comb_intersection_length_list.append(len(set(cat_lists[0]).intersection(*cat_lists)))

# removendo combinacoes de categorias com intersecoes de tamanho 0
comb_list_list = np.array(comb_list_list)
comb_intersection_length_list = np.array(comb_intersection_length_list)
comb_list_list = comb_list_list[comb_intersection_length_list != 0]
comb_intersection_length_list = comb_intersection_length_list[comb_intersection_length_list != 0]

print(sorted(comb_intersection_length_list,reverse=True))
# criando uma serie de membros que indica o tamanho da intersecao entre os diferentes conjuntos
mem_series = upsetplot.from_memberships(comb_list_list, data=comb_intersection_length_list)

# plot
upsetplot.plot(mem_series, orientation='horizontal',sort_by='cardinality', show_counts=False, totals_plot_elements=False, show_percentages=False)

plt.grid(False)

# save
plt.savefig("upset_plot.png", bbox_inches='tight')
#plt.show()
