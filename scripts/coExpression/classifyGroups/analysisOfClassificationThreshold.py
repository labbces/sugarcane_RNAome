#!/usr/bin/env python3

from collections import defaultdict
import matplotlib.pyplot as plt

def classify_group(classifications, threshold):
    total_genes = len(classifications)
    protein_coding_genes = classifications.count('protein-coding')
    non_coding_genes = classifications.count('non-coding')
    protein_and_non_coding_genes = classifications.count('protein and non-coding')
    
    protein_coding_percentage = (protein_coding_genes / total_genes) * 100
    non_coding_percentage = (non_coding_genes / total_genes) * 100
    protein_and_non_coding_percentage = (protein_and_non_coding_genes / total_genes) * 100

    if protein_coding_percentage >= threshold:
        return 'protein-coding'
    elif non_coding_percentage >= threshold:
        return 'non-coding'
    elif protein_and_non_coding_percentage >= threshold:
        return 'protein and non-coding'
    else:
        return 'unknown'

group_classifications = defaultdict(list)

with open('panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv', 'r') as file:
    for line in file:
        line = line.strip().split('\t')
        group_id = line[1]
        classification = line[3]

        group_classifications[group_id].append(classification)

thresholds = [50, 60, 70, 80, 90]

for threshold in thresholds:
    classification_counts = defaultdict(int)

    for group_id, classifications in group_classifications.items():
        group_classification = classify_group(classifications, threshold)
        classification_counts[group_classification] += 1

    output_filename_txt = f'classificationByPanRNAomeGroup_{threshold}.txt'
    with open(output_filename_txt, 'w') as output_file:
        output_file.write("Classification\tCount\n")
        for classification, count in classification_counts.items():
            output_file.write(f"{classification}\t{count}\n")

    labels, values = zip(*classification_counts.items())
    fig, ax = plt.subplots()
    ax.bar(labels, values, color=['blue', 'green', 'orange', 'gray'])
    ax.set_ylabel('Count')
    ax.set_title(f'Distribution of Classifications by PanRNAome Group ({threshold}% threshold)')
    plt.tight_layout()

    output_filename_png = f'classificationByPanRNAomeGroup_{threshold}.png'
    plt.savefig(output_filename_png)
