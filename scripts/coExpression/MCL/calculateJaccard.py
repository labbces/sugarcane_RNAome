#!/usr/bin/env python3

import os
import argparse
from itertools import combinations

parser = argparse.ArgumentParser(description="Calculate similarities between clusters from MCL")
parser.add_argument("-d", dest="directory", type=str, metavar="clusters_directory", help="Directory containing clusters from MCL", required=True)
parser.add_argument("-o", dest="output", type=str, metavar="output_clusteringSimilarities.tsv", help="output name", required=True)
parser.add_argument("-t", dest="threshold", type=float, metavar="threshold < int >", help="0.5", required=True)

args = parser.parse_args()
directory = args.directory
output = args.output + "_clusteringSimilarities.tsv"
threshold = args.threshold

# load clusters from file
def load_clusters(filename):
    with open(filename, 'r') as file:
        #clusters = [line.strip().split('\t') for line in file.readlines()]
        clusters = []
        for line in file:
            if line.startswith("clusters="):
                break
            else:
                clusters.append(line.strip().split('\t'))
    return clusters

# calculate Jaccard coefficient
def jaccard(cluster1, cluster2):
    intersection = len(set(cluster1).intersection(cluster2))
    union = len(set(cluster1).union(cluster2))
    return intersection / union if union != 0 else 0

# calculate Overlap coefficient
def overlap(cluster1, cluster2):
    intersection = len(set(cluster1).intersection(cluster2))
    min_size = min(len(cluster1), len(cluster2))
    return intersection / min_size if min_size != 0 else 0

# Load clusters from each set
cluster_sets = []
cluster_filenames = []

for filename in os.listdir(directory):
    clusters = load_clusters(os.path.join(directory, filename))
    #print(f'clusters: {clusters}')
    cluster_sets.append(clusters)
    cluster_filenames.append(filename)  # Remove extension from filename

total_calculations = 0
total_combinations = sum(len(cluster_sets[i]) * len(cluster_sets[j]) for i, j in combinations(range(len(cluster_sets)), 2))

print(f'Total clusters combinations: {total_combinations}')
print(f'Saving clusters with similarity > {threshold}\n')

# Write the results to output file
with open(output, "w") as output_file:
    output_file.write("cluster_1\tcluster_2\tJaccard\tOverlap\n")
    for i, (clusters1, filename1) in enumerate(zip(cluster_sets, cluster_filenames), start=1):
        for j, (clusters2, filename2) in enumerate(zip(cluster_sets[i:], cluster_filenames[i:]), start=i+1):
            for idx1, cluster1 in enumerate(clusters1, start=1):
                for idx2, cluster2 in enumerate(clusters2, start=1):
                    
                    jaccard_coefficient = jaccard(cluster1, cluster2)
                    overlap_coefficient = overlap(cluster1, cluster2)
                    
                    total_calculations += 1
                    
                    if jaccard_coefficient > threshold or overlap_coefficient > threshold:
                        output_file.write(f"{filename1}_{idx1}\t{filename2}_{idx2}\t{jaccard_coefficient}\t{overlap_coefficient}\n")

                    progress_percentage = total_calculations / total_combinations * 100
                    print(f"Similarities coefficient calculated for clusters: {total_calculations}/{total_combinations} - {progress_percentage:.2f}%", end="\r")
                    
print(f'\nResults saved to {output}')