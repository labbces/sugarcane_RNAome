#!/usr/bin/env python3

import pandas as pd

classification_file = "panTranscriptome_panRNAomeClassificationTable_hyphen.tsv"
ncRNA_file = "putative_ncRNA_consensus.txt"
proteinCoding_file = "proteinCoding_RNA.txt"

with open(ncRNA_file, 'r') as f:
    ncRNA_transcripts = set(f.read().splitlines())

with open(proteinCoding_file, 'r') as f:
    proteinCoding_transcripts = set(f.read().splitlines())

df = pd.read_csv(classification_file, sep='\t', header=None)
# Soft-core	OG0000000	B1_k25_TRINITY_DN12555_c1_g1_i10
df.columns = ['Category', 'Gene', 'Transcript']

gene_classification = {}

for index, row in df.iterrows():
    gene = row['Gene']
    transcript = row['Transcript']
    
    if transcript in proteinCoding_transcripts and transcript in ncRNA_transcripts:
        classification = 'protein and non-coding'
    elif transcript in proteinCoding_transcripts:
        classification = 'protein-coding'
    elif transcript in ncRNA_transcripts:
        classification = 'non-coding'
    #else:
    #    classification = 'unknown'
    
    if gene not in gene_classification:
        gene_classification[gene] = classification
    else:
        if gene_classification[gene] != classification:
            gene_classification[gene] = 'protein and non-coding'

df['Gene_Classification'] = df['Gene'].map(gene_classification)

output_file = "panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv"
df.to_csv(output_file, sep='\t', index=False, header=False)
