#!/usr/bin/env python3

import os
import argparse

parser = argparse.ArgumentParser(description="Extract genes from consistent clusters - cliques")
parser.add_argument("-d", dest="clusters_dir", type=str, metavar="clusters_directory", help="Directory containing clusters from MCL", required=True)
parser.add_argument("-o", dest="output_dir", type=str, metavar="output_directory", help="cliques", required=True)
parser.add_argument("-c", dest="cliques_file", type=str, metavar="cliques_file", help="cliques_jaccard.txt", required=True)

args = parser.parse_args()
clusters_dir = args.clusters_dir 
cliques_file = args.cliques_file
output_dir = args.output_dir

#cliques_file = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/correr/MCL/CNC/clusterSimilarity/cliques_jaccard.txt"
#clusters_dir = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/correr/MCL/CNC/"
#output_dir = "/home/felipe/Documents/sugarcane_RNAome/scripts/coExpression/fiberAndSugar/correr/MCL/CNC/Cliques/"

cliques = {}

with open(cliques_file, 'r') as cliques_file:
    for idx, line in enumerate(cliques_file, start=1):
        clique_name = f"clique_{idx}"
        cluster_members = line.strip().split()
        cliques[clique_name] = cluster_members
        
#print(cliques)

# make sure outdir exists
os.makedirs(output_dir, exist_ok=True)

# Test - saving only one clique file
#clique_to_save = "clique_5"

for clique, clusters in cliques.items():
    
    # just for testing - only one clique
    #if clique == clique_to_save:
    
    clique_genes = set()
    gene_clusters = {}

    # for each cluster file in this clique
    for cluster_info in clusters:
        # extract cluster name and line
        cluster_file, line_number = cluster_info.split("_")

        cluster_path = os.path.join(clusters_dir, cluster_file)

        if os.path.exists(cluster_path):
            # Read each cluster file and add genes from specific line into clique_genes set
            with open(cluster_path, 'r') as file:
                for idx, line in enumerate(file, start=1):
                    if idx == int(line_number):
                        genes = line.strip().split()
                        clique_genes.update(genes)
                        
                        # Also update cluster dict with clique genes
                        for gene in genes:
                            if gene in gene_clusters:
                                gene_clusters[gene].append((cluster_file, line_number))
                            else:
                                gene_clusters[gene] = [(cluster_file, line_number)]
                        break
        else:
            print(f"Cluster file {cluster_file} not found.")

    # Save genes and clusters for each clique
    output_file = os.path.join(output_dir, f"{clique}.txt")
    
    with open(output_file, 'w') as output:
        for gene in clique_genes:
            # if the gene is present in all clusters: "intersection"
            # else: "disjoint"
            
            classification = "intersection" if len(gene_clusters[gene]) == len(clusters) else "disjoint"
            
            member_of = ", ".join([f"{cluster}_{line}" for cluster, line in gene_clusters[gene]])
            
            output.write(f"{gene}\t{classification}\t{member_of}\n")
            
    print(f"Genes from clique {clique} saved in {output_file}!")