#!/usr/bin/env python3

# Classify individual group
def classify_group(group_data):
    classifications = set()
    for item in group_data:
        if "protein-coding" in item:
            classifications.add("protein-coding")
        elif "non-coding" in item:
            classifications.add("non-coding")
        else:
            classifications.add("protein and non-coding")
    if len(classifications) == 1:
        return classifications.pop()
    elif "protein-coding" in classifications and "non-coding" in classifications:
        return "protein and non-coding"
    else:
        return "unknown"  

group_results = {}

infile = "panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv"

with open(infile, 'r') as file:
    for line in file:
        data = line.strip().split('\t')
        group_id = data[1]  # "Group_ID"
        
        if group_id not in group_results:
            group_results[group_id] = {"Group": data[0], "Classification": [data[3]]}
        else:
            group_results[group_id]["Classification"].append(data[3])

with open('CNC_groups.txt', 'w') as output_file:
    for group_id, result in group_results.items():
        group_classification = classify_group(result["Classification"])
        output_file.write(f"{result['Group']}\t{group_id}\t{group_classification}\n")

print("Resultados foram escritos em 'CNC_groups.txt'.")
