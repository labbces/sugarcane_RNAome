import pandas as pd
import matplotlib.pyplot as plt

# Ler o arquivo com os dados dos graus dos nós
file_path = "Perlo2022_counts_filters_VST_CNC_CV_above0.6_mcl_degree_sorted.tsv"
data = pd.read_csv(file_path, delimiter="\t", header=None, names=["Node", "Degree"])

# Plotar a distribuição de links para os nós com diferentes números de bins
bin_sizes = [10, 100, 500, 1000]
plt.figure(figsize=(15, 10))

for i, bins in enumerate(bin_sizes, 1):
    plt.subplot(2, 2, i)
    plt.hist(data["Degree"], bins=bins, color='black', edgecolor='black')
    plt.title(f"Degree distribution - {bins} Bins")
    plt.xlabel("Degree")
    plt.ylabel("Nodes")
    plt.grid(True)

plt.tight_layout()
plt.show()
