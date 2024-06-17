import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

file_paths_titles = [
    ("Perlo2022_counts_filters_VST_CNC_CV_above0.6_mcl_degree_sorted_updated.tsv", "CDF of degree distribution - Perlo2022 - CV > 0.6"),
    ("Correr2020_counts_filters_VST_CNC_CV_above_1.5_mcl_combined_degree_sorted_updated.tsv", "CDF of degree distribution - Correr2020 - CV > 1.5"),
    ("Correr2020_counts_filters_VST_CNC_CV_above2_mcl_degree_sorted_updated.tsv", "CDF of degree distribution - Correr2020 - CV > 2.0"),
    ("Hoang2017_counts_filters_VST_CNC_CV_above_1_mcl_degree_sorted_updated.tsv", "CDF of degree distribution - Hoang2017 - CV > 1.0"),
    ("Hoang2017_counts_filters_VST_CNC_CV_above1.2_mcl_degree_sorted_updated.tsv", "CDF of degree distribution - Hoang2017 - CV > 1.2")
]

colors = {'lncRNA': 'green', 'ncRNA': 'blue', 'protein-coding': 'red'}

def compute_complementary_cdf(data):
    counts, bin_edges = np.histogram(data, bins=len(data), density=True)
    cdf = np.cumsum(counts) / np.sum(counts)
    return bin_edges[1:], 1 - cdf  # 1 - CDF for complementary CDF

''' Notes for the future:
The CDF represents the cumulative probability of the degree distribution. 
By definition, it starts from 0 and goes up to 1. When we plot it on a logarithmic scale, the values are transformed. 
Since the values are between 0 and 1, on a logarithmic scale, they start at 10^-1 (which is 0.1) and go up to 10^0 (which is 1).

The reference article shows the complementary CDF P(X>k), which represents the probability that a node's degree is greater than k. 
This is a decreasing function because as k increases, fewer nodes have degrees greater than k.
'''

for file_path, title in file_paths_titles:

    data = pd.read_csv(file_path, delimiter="\t", header=None, names=["Node", "Degree", "panRNAome classification", "Transcript", "Gene function", "Transcript size", "Transcript function"])
    
    plt.figure(figsize=(10, 6))

    for transcript_function in colors.keys():
        subset = data[data["Transcript function"] == transcript_function]
        degrees = subset["Degree"]
        x, y = compute_complementary_cdf(degrees)
        plt.plot(x, y, 'o', label=transcript_function, color=colors[transcript_function])
    
    plt.xscale('log')
    plt.yscale('log')
    plt.title(f"{title}")
    plt.xlabel("Degree(k)")
    plt.ylabel("P(K>k)")
    plt.grid(True, which="both", ls="--")
    plt.legend(title='Transcript function')
    #plt.show()

    output_file = os.path.splitext(file_path)[0] + "_cdf_degree_distribution.png"
    
    plt.savefig(output_file)
    plt.close()