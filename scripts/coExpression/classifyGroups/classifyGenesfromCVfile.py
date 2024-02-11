#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description="Script to classify function and panRNAome groups of filtered expression matrix and plot histogram.")
parser.add_argument("-c", dest="classification_file", metavar="panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv", required=True, help="Path to the classification file")
parser.add_argument("-cv", dest="CV_file", metavar="top20CV.txt", required=True, help="Path to the CV file")
parser.add_argument("-o", dest="output_file", metavar="top20CV_classified.txt", required=True, help="Path to the output classified file")
parser.add_argument("-opng", dest="output_filename_png", metavar="top % CV file classified PNG", required=False, help="Path to the output PNG file for the histogram")

args = parser.parse_args()
classification_file = args.classification_file
cv_file = args.CV_file
out = args.output_file
output_filename_png = args.output_filename_png

gene_functions = {}
gene_classifications = {}

with open(classification_file, 'r') as file:
    for line in file:
        fields = line.strip().split('\t')
	#fields = line.strip().split() # if sep is space not tab
        classification = fields[0]
        gene = fields[1]
        transcript = fields[2]
        function = fields[3]
        if gene not in gene_functions:
            gene_functions[gene] = set()
        if gene not in gene_classifications:
            gene_classifications[gene] = set()
        gene_functions[gene].add(function)
        gene_classifications[gene].add(classification)

def get_predominant_function(gene):
    functions = gene_functions.get(gene, set())
    if "protein and non-coding" in functions:
    #if "protein" in functions: # if sep is space not tab
        return "protein-coding"
    elif "protein-coding" in functions:
        return "protein-coding"
    elif "non-coding" in functions:
        return "non-coding"

def get_predominant_classification(gene):
    classifications = gene_classifications.get(gene, set())
    if "Accessory" in classifications:
        return "Accessory"
    elif "Soft-core" in classifications:
        return "Soft-core"
    elif "Hard-core" in classifications:
        return "Hard-core"
    elif "Exclusive" in classifications:
        return "Exclusive"

with open(cv_file, 'r') as file:
    with open(out, 'w') as output_file:
        for line in file:
            gene, cv = line.strip().split('\t')
            function = get_predominant_function(gene)
            classification = get_predominant_classification(gene)
            output_file.write(f"{gene}\t{cv}\t{classification}\t{function}\n")
