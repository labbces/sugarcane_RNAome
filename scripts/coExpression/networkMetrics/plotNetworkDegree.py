import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os 

file_paths_titles = [
    ("Perlo2022_counts_filters_VST_CNC_CV_above0.6_mcl_degree_sorted_updated.tsv", "Degree distribution - Perlo2022 - CV > 0.6 - {bins} Bins"),
    ("Correr2020_counts_filters_VST_CNC_CV_above_1.5_mcl_combined_degree_sorted_updated.tsv", "Degree distribution - Correr2020 - CV > 1.5 - {bins} Bins"),
    ("Correr2020_counts_filters_VST_CNC_CV_above2_mcl_degree_sorted_updated.tsv", "Degree distribution - Correr2020 - CV > 2.0 - {bins} Bins"),
    ("Hoang2017_counts_filters_VST_CNC_CV_above_1_mcl_degree_sorted_updated.tsv", "Degree distribution - Hoang2017 - CV > 1.0 - {bins} Bins"),
    ("Hoang2017_counts_filters_VST_CNC_CV_above1.2_mcl_degree_sorted_updated.tsv", "Degree distribution - Hoang2017 - CV > 1.2 - {bins} Bins")
]

colors = {'lncRNA': 'green', 'ncRNA': 'blue', 'protein-coding': 'red'}

bin_sizes = [10, 100, 500, 1000]

for file_path, title_template in file_paths_titles:
    data = pd.read_csv(file_path, delimiter="\t", header=None, names=["Gene", "Degree", "panRNAome classification", "Transcript", "Gene function", "Transcript size", "Transcript function"])
    
    plt.figure(figsize=(15, 10))

    for i, bins in enumerate(bin_sizes, 1):
        plt.subplot(2, 2, i)
        
        sns.histplot(data, x='Degree', hue='Transcript function', bins=bins, palette=colors, multiple='stack', edgecolor='black')
        
        plt.title(title_template.format(bins=bins))
        plt.xlabel("Degree")
        plt.ylabel("Nodes")
        plt.grid(True)

    plt.tight_layout()
    #plt.show()

    output_file = os.path.splitext(file_path)[0] + "_degree_distribution.png"
    
    plt.savefig(output_file)
    plt.close()
