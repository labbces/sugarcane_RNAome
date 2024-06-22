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

# transcript function
#colors = {'lncRNA': 'green', 'ncRNA': 'blue', 'protein-coding': 'red'}

# panRNAome category 
colors = {'Exclusive': 'green', 'Accessory': 'blue', 'Soft-core': 'red', 'Hard-core': 'yellow'}

bin_sizes = [10, 100, 500, 1000]

for file_path, title_template in file_paths_titles:
    data = pd.read_csv(file_path, delimiter="\t", header=None, names=["Gene", "Degree", "panRNAome classification", "Transcript", "Gene function", "Transcript size", "Transcript function"])
    
    # filter transcript by function
    lncrnas = data.loc[data['Transcript function']=='lncRNA']
    ncrnas = data.loc[data['Transcript function']=='ncRNA']
    protein_coding = data.loc[data['Transcript function']=='protein-coding']
    
    plt.figure(figsize=(15, 10))

    for i, bins in enumerate(bin_sizes, 1):
        plt.subplot(2, 2, i)
        
        # plot degree by transcript function
        #sns.histplot(data, x='Degree', hue='Transcript function', bins=bins, palette=colors, multiple='stack', edgecolor='black')
        
        # filter and plot degree by transcript function -> focusing on panRNAome categories
        sns.histplot(lncrnas, x='Degree', hue='panRNAome classification', bins=bins, palette=colors, multiple='stack', edgecolor='black')
        #sns.histplot(ncrnas, x='Degree', hue='panRNAome classification', bins=bins, palette=colors, multiple='stack', edgecolor='black')
        #sns.histplot(protein_coding, x='Degree', hue='panRNAome classification', bins=bins, palette=colors, multiple='stack', edgecolor='black')

        plt.title(title_template.format(bins=bins))
        plt.xlabel("Degree")
        plt.ylabel("Nodes")
        plt.grid(True)

    plt.tight_layout()
    #plt.show()

    # change outname for each version (transcript function and panRNAome category)
    #output_file = os.path.splitext(file_path)[0] + "_degree_distribution.png"
    output_file = os.path.splitext(file_path)[0] + "_degree_distribution_lncRNA_category.png"

    plt.savefig(output_file)
    plt.close()
