#!/usr/bin/env python3

from collections import defaultdict

# Leitura do arquivo de entrada
input_file = "DB_clust.tsv"
output_file = "DB_clust_Orthogroups.tsv"

orthogroups = defaultdict(lambda: defaultdict(list))
genotypes = set()

# Leitura do arquivo e processamento dos dados
with open(input_file, "r") as file:
    for line in file:
        columns = line.strip().split()
        if len(columns) == 2:  # Linhas com dois elementos representam relações
            orthogroup, element = columns
            genotype = element.split('_')[0]  # Obtém o genótipo
            orthogroups[orthogroup][genotype].append(element)
            genotypes.add(genotype)

# Ordena os genótipos alfabeticamente
sorted_genotypes = sorted(genotypes)

# Criação do novo formato e impressão
output_lines = ["Orthogroup\t" + "\t".join(sorted_genotypes)]

for idx, (orthogroup, genotypes_dict) in enumerate(orthogroups.items(), start=1):
    genotype_elements = []
    for genotype in sorted_genotypes:
        elements = genotypes_dict.get(genotype, [])
        elements_str = ",".join(elements)
        genotype_elements.append(elements_str)
    output_lines.append("OG{}\t{}".format(idx, "\t".join(genotype_elements)))

with open(output_file, "w") as file:
    file.write("\n".join(output_lines))

print("Arquivo de saída criado com sucesso:", output_file)

