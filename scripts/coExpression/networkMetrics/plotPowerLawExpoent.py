import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import powerlaw

file_path = "Perlo2022_counts_filters_VST_CNC_CV_above0.6_mcl_degree_sorted.tsv"
data = pd.read_csv(file_path, delimiter="\t", header=None, names=["Node", "Degree"])

# calculate degree distribution
degree_counts = data["Degree"].value_counts().sort_index()

# plot degree distribution
plt.figure(figsize=(10, 6))
plt.loglog(degree_counts.index, degree_counts.values, 'bo', markersize=5)
plt.title("Degree distribution (Log-Log)")
plt.xlabel("Degree")
plt.ylabel("Nodes")
plt.grid(True)
plt.show()

# adjust power law expoent
degree_values = degree_counts.index
frequency_values = degree_counts.values

# remove zeros (dealing with log transformation)
degree_values = degree_values[frequency_values > 0]
frequency_values = frequency_values[frequency_values > 0]

# adjust line to log-log
log_degree = np.log(degree_values)
log_frequency = np.log(frequency_values)

slope, intercept = np.polyfit(log_degree, log_frequency, 1)
gamma = -slope

plt.figure(figsize=(10, 6))
plt.loglog(degree_values, frequency_values, 'bo', markersize=5, label="Nodes")
plt.plot(degree_values, np.exp(intercept) * degree_values**slope, 'r', linewidth=2, label=f"P(k) = {np.exp(intercept):.2f} * k^{{-{gamma:.2f}}}")
plt.title("Degree distribution (Log-Log)")
plt.xlabel("Degree")
plt.ylabel("Nodes")
plt.grid(True)
plt.legend()
plt.show()

# power law exponent
print(f"power law exponent: {gamma:.2f}")