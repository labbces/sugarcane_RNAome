import pandas as pd
import matplotlib.pyplot as plt

columns = ["clusters", "max", "ctr", "avg", "min", "DGI", "TWI", "TWL", "sgl", "qrt", "efficiency", "massfrac", "areafrac", "target-name", "source-name"]

file = pd.read_csv("../fiberAndSugar/hoang/MCL/CNC/CV1.0/clminfo.csv", sep=",", usecols=columns, index_col=False)
file["inflation"] = file["source-name"].str.extract("(\d+\.\d+)")

print(file["inflation"])
print(file["source-name"])

fig, axs = plt.subplots(1, 3, figsize=(12, 5)) # create figure with subplots

# Plot 1 - efficiency vs inflation
axs[0].plot(file["efficiency"], file["inflation"], marker="o",  label="efficiency", color="black")
axs[0].set_xlabel("Efficiency")
axs[0].set_ylabel("Inflation")
axs[0].set_title("Clustering Efficiency vs Inflation")
axs[0].legend(frameon=False)
axs[0].grid(False)

# Plot 3 - cluster size, max cluster size vs inflation
axs[1].plot(file["clusters"], file["inflation"], marker="o", label="clusters", color="black")
axs[1].plot(file["max"], file["inflation"], marker="^", label="max size", color="darkgrey")
axs[1].set_xlabel("Clusters")
#axs[1].set_ylabel("Inflation")
axs[1].set_title("Cluster size vs Inflation")
axs[1].legend(frameon=False)
axs[1].grid(False)

# Plot 4 - mass fraction - area fraction vs inflation
axs[2].plot(file["massfrac"], file["inflation"], marker="o", label="mass fraction", color="black")
axs[2].plot(file["areafrac"], file["inflation"], marker="^", label="area fraction", color="darkgrey")
axs[2].set_xlabel("Fraction")
#axs[2].set_ylabel("Inflation")
axs[2].set_title("Fraction vs Inflation")
axs[2].legend(frameon=False)
axs[2].grid(False)

# Hide y-axis for the second and third subplots
#axs[1].get_yaxis().set_visible(False)
#axs[2].get_yaxis().set_visible(False)

plt.tight_layout() # adjust layout

plt.savefig("clusteringEfficiency.png")
plt.show()