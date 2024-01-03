#!/usr/bin/env python3

with open('proteinCoding_RNA.txt', 'r') as protein_file:
    protein_coding_genes = set(line.strip() for line in protein_file)

with open('putative_ncRNA_consensus.txt', 'r') as non_coding_file:
    non_coding_genes = set(line.strip() for line in non_coding_file)

common_genes = protein_coding_genes.intersection(non_coding_genes)

coding_non_coding = 0
protein_coding = 0
non_coding = 0
unknown = 0

classified_genes = {}

with open('panTranscriptome_panRNAomeClassificationTable_hyphen.tsv', 'r') as input_file, open('panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv', 'w') as output_file:
    for line in input_file:
        gene = line.split('\t')[2].strip()

        if gene in common_genes:
            classification = 'protein and non-coding'
            coding_non_coding += 1 

        elif gene in protein_coding_genes:
            classification = 'protein-coding'
            protein_coding +=1        

        elif gene in non_coding_genes:
            classification = 'non-coding'
            non_coding += 1

        else:
            classification = 'unknown'
            unknown += 1
        
        if gene not in classified_genes:
            classified_genes[gene] = {'classification': classification, 'count': 0}
      
        classified_genes[gene]['count'] += 1
        output_file.write(line.strip() + '\t' + classification + '\n')

def print_metric(classification, label):
    unique_genes = set(gene for gene, info in classified_genes.items() if info['classification'] == classification)
    print(f"{label}: {len(unique_genes)}")

print("### FINISHED GROUPS CLASSIFICATION")
print("\n")
print_metric('protein and non-coding', 'protein and non-coding sequences')
print_metric('protein-coding', 'protein-coding sequences')
print_metric('non-coding', 'non-coding sequences')
print_metric('unknown', 'unknown sequences')

print('\n')
print("Total classified genes:", len(classified_genes))
print("Total lines processed (classified_genes['count']):", sum(info['count'] for info in classified_genes.values()))

print("\n")

total_duplicates = sum(info['count'] -1 for info in classified_genes.values() if info['count'] > 1)
print(f"Total duplicates: {total_duplicates}")
