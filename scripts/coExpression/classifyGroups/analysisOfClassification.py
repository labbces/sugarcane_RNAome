#!/usr/bin/env python3

import matplotlib.pyplot as plt
from collections import defaultdict

group_classifications = defaultdict(set)

with open('panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv', 'r') as file:
    for line in file:
        line = line.strip().split('\t')
        group_id = line[1]
        classification = line[3]

        group_classifications[group_id].add(classification)

classification_counts = defaultdict(int)

for group_id, classifications in group_classifications.items():
    if 'protein and non-coding' in classifications:
        classification_counts['protein and non-coding'] += 1
    elif 'protein-coding' in classifications:
        classification_counts['protein-coding'] += 1
    elif 'non-coding' in classifications:
        classification_counts['non-coding'] += 1
    else:
        classification_counts['unknown'] += 1

with open('classificationByPanRNAomeGroup.txt', 'w') as output_file:
    output_file.write("Classification\tCount\n")
    for classification, count in classification_counts.items():
        output_file.write(f"{classification}\t{count}\n")

# Plot histogram
labels, values = zip(*classification_counts.items())
fig, ax = plt.subplots()
ax.bar(labels, values, color=['blue', 'green', 'orange', 'gray'])
ax.set_ylabel('Count')
ax.set_title('Distribution of Classifications by panRNAome Groups')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

# Salvando o plot em um arquivo
plt.savefig('classificationByPanRNAomeGroup.png')
