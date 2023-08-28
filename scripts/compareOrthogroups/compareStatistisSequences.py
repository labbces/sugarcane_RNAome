import matplotlib.pyplot as plt

def read_metrics_file(filename):
    metrics = {}
    with open(filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            key, value = line.strip().split(":")
            metrics[key] = float(value)
    return metrics

protein_metrics = read_metrics_file('proteinCodingPanTranscriptomeStatisticsOrthogroups_0.8.tsv')
noncoding_metrics = read_metrics_file('nonCodingPanTranscriptomeStatisticsOrthogroups_0.8.tsv')
transcript_metrics = read_metrics_file('panTranscriptomeStatisticsOrthogroups_0.8.tsv')

metric_names = [
    #"Number of groups present in all (48) species",
    "Number of proteins present in core groups",
    #"Number of groups present in 80% of the species (38.400000000000006)",
    "Number of proteins present in soft-core groups",
    #"Number of accesory groups",
    "Number of proteins present in accessory groups",
    #"Number of exclusive groups"
    "Number of proteins present in exclusive groups"
]

# Prepare data for plotting
data = {
    'protein-coding': [protein_metrics[metric] for metric in metric_names],
    'non-coding': [noncoding_metrics[metric] for metric in metric_names],
    'pan-transcriptoma': [transcript_metrics[metric] for metric in metric_names]
}

# Create a bar plot
fig, ax = plt.subplots(figsize=(10, 6))
width = 0.2
x = range(len(metric_names))

for idx, (label, values) in enumerate(data.items()):
    ax.bar(
        [pos + width * idx for pos in x],
        values,
        width=width,
        label=label
    )

ax.set_xticks([pos + width for pos in x])
ax.set_xticklabels(metric_names, rotation=45, ha="right")
ax.set_ylabel('Sequences (Log Scale)')
ax.set_title('Comparison of Sequences between soft-core values (80%)')
ax.legend()

# Set the y-axis to log scale and adjust y-axis limits for better visualization
ax.set_yscale('log')
ax.set_ylim(bottom=0.1)  # Adjust this value as needed

plt.savefig("comparisonSoft-CoreSequences_0.8.png", bbox_inches="tight")
plt.tight_layout()
plt.show()