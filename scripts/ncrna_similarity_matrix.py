#!/usr/bin/env python

# Modificar DB_clust.tsv
# awk -F '\t' 'BEGIN {OFS = FS} {split($2, arr, "_"); $2 = arr[1]; print}' DB_clust.tsv >> DB_clust_genotypeName.tsv 

# Usage:
#./ncrna_similarity_matrix.py > DB_clust_groups.tsv

from collections import defaultdict

# Leitura do arquivo de similaridades
similarities = []
with open('DB_clust_genotypeName.tsv', 'r') as file:
#with open('teste.tsv', 'r') as file:
    for line in file:
        gene1, gene2 = line.strip().split('\t')
        similarities.append((gene1, gene2))

# Construção da matriz de grupos
groups = defaultdict(dict)
group_count = 1  # Contador para gerar os nomes dos grupos
for gene1, gene2 in similarities:
    if gene1 not in groups:
        group_name = 'OG' + str(group_count)
        group_count += 1
        groups[gene1]['group_name'] = group_name
    genotype = gene2
    if genotype not in groups[gene1]:
        groups[gene1][genotype] = 1
    else:
        groups[gene1][genotype] += 1

# Obter todos os genótipos
genotypes = set()
for group in groups.values():
    genotypes.update(group.keys())
genotypes = sorted(genotypes)

# Impressão da matriz
print("Orthogroup", end='\t')
print('\t'.join(genotypes))

for gene1, group in groups.items():
    group_name = group['group_name']
    print(group_name, end='\t')
    counts = [str(group.get(genotype, 0)) for genotype in genotypes]
    print('\t'.join(counts))
