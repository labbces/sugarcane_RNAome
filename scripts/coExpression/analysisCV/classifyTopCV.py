#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Script to classify panRNAome groups of top % CV and plot a histogram.")
parser.add_argument("-c", dest="classification_file", metavar="panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv", required=True, help="Path to the classification file")
parser.add_argument("-g", dest="group_file", metavar="top20CV.txt", required=True, help="Path to the group file")
parser.add_argument("-o", dest="output_file", metavar="top20CV_classified.txt", required=True, help="Path to the output classified file")
parser.add_argument("-opng", dest="output_filename_png", metavar="top % CV file classified PNG", required=True, help="Path to the output PNG file for the histogram")
args = parser.parse_args()

classifications = {}
with open(args.classification_file, 'r') as f:
    for line in f:
        fields = line.strip().split('\t')
        group_id = fields[1]
        classification = fields[3]
        classifications[group_id] = classification

classifications_counts = {}
with open(args.group_file, 'r') as f, open(args.output_file, 'w') as out_file:
    for line in f:
        fields = line.strip().split('\t')  # fields[0] = group // fields[1] = cv
        if len(fields) > 1:  # Ignore lines with less than 2 columns
            group_id = fields[0]
            if group_id in classifications:
                classification = classifications[group_id]
                if classification in classifications_counts:
                    classifications_counts[classification] += 1
                else:
                    classifications_counts[classification] = 1
                out_file.write(f"{group_id}\t{classification}\t{fields[1]}\n")

print("Classification of CV file done.", "\n")

print("Plotting histogram ...", "\n")

labels, values = zip(*classifications_counts.items())
colors = ['gold', 'lightcoral', 'lightskyblue', 'lightgreen'] 

plt.figure(figsize=(10, 6))
plt.bar(labels, values, color=colors)
plt.xlabel("Classification")
plt.ylabel("Count")
plt.title("Distribution of Classifications by PanRNAome Group (top 20% CV)")

plt.savefig(args.output_filename_png)
print(f"Histogram saved as {args.output_filename_png}")

