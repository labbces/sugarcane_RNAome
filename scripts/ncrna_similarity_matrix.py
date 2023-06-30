from collections import defaultdict

# Função para ler o arquivo TSV de similaridades e retornar uma lista de pares de genes
def ler_lista_similaridades(arquivo):
    similaridades = []
    with open(arquivo, "r") as file:
        for line in file:
            gene1, gene2 = line.strip().split("\t")
            similaridades.append([gene1, gene2])
    return similaridades

# Função para transformar a lista de similaridades em uma matriz de orthogroups
def gerar_matriz_orthogroups(similaridades):
    gene_to_orthogroup = {}
    orthogroup_to_genes = defaultdict(list)

    for genes in similaridades:
        gene1, gene2 = genes
        if gene1 not in gene_to_orthogroup and gene2 not in gene_to_orthogroup:
            orthogroup = "OG" + str(len(gene_to_orthogroup) + 1)
            gene_to_orthogroup[gene1] = orthogroup
            gene_to_orthogroup[gene2] = orthogroup
            orthogroup_to_genes[orthogroup].extend([gene1, gene2])
        elif gene1 in gene_to_orthogroup and gene2 not in gene_to_orthogroup:
            orthogroup = gene_to_orthogroup[gene1]
            gene_to_orthogroup[gene2] = orthogroup
            orthogroup_to_genes[orthogroup].append(gene2)
        elif gene2 in gene_to_orthogroup and gene1 not in gene_to_orthogroup:
            orthogroup = gene_to_orthogroup[gene2]
            gene_to_orthogroup[gene1] = orthogroup
            orthogroup_to_genes[orthogroup].append(gene1)
        else:
            orthogroup1 = gene_to_orthogroup[gene1]
            orthogroup2 = gene_to_orthogroup[gene2]
            if orthogroup1 != orthogroup2:
                orthogroup_to_genes[orthogroup1].extend(orthogroup_to_genes[orthogroup2])
                for gene in orthogroup_to_genes[orthogroup2]:
                    gene_to_orthogroup[gene] = orthogroup1
                del orthogroup_to_genes[orthogroup2]

    orthogroups = []
    for orthogroup, genes in orthogroup_to_genes.items():
        orthogroups.append([orthogroup] + sorted(genes))

    return sorted(orthogroups)

# Função para imprimir a matriz de orthogroups
def imprimir_matriz_orthogroups(orthogroups):
    for orthogroup in orthogroups:
        print("\t".join(orthogroup))

# Caminho para o arquivo TSV de similaridades
arquivo_similaridades = "caminho/para/o/arquivo.tsv"

# Ler a lista de similaridades a partir do arquivo TSV
similaridades = ler_lista_similaridades(arquivo_similaridades)

# Gerar a matriz de orthogroups
orthogroups = gerar_matriz_orthogroups(similaridades)

# Imprimir a matriz de orthogroups
imprimir_matriz_orthogroups(orthogroups)
