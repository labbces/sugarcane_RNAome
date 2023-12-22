#!/usr/bin/env python3

def load_table(file_path):
    data = set()
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            if len(columns) >= 3:
                data.add(columns[2])
    return data

def load_quant(file_path):
    data = set()
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            data.add(columns[0])
    return data

def find_exact_matches(table1, table2):
    exact_matches = table1.intersection(table2)
    return exact_matches

file_path1 = "SRR5258946_quant.sf"
file_path2 = "panTranscriptome_panRNAomeClassificationTable.tsv"

# load tables
table1 = load_quant(file_path1)
table2 = load_table(file_path2)

exact_matches = find_exact_matches(table1, table2)
  
output_file_path = "panTranscriptome_panRNAomeClassificationTable_smallData.tsv"

with open(file_path2, 'r') as file2,  open(output_file_path, 'w') as output_file:
    output_file.write("sequence matches in SRR5258946_quant.sf and panTranscriptome_panRNAomeClassificationTable.tsv:\n")
    for line in file2:
        columns = line.strip().split('\t')
        if columns[2] in exact_matches:
            output_file.write(columns[0] + '\t' + columns[1] + '\t' + columns[2] + '\n')

print(f"Exact matches saved in: {output_file_path}")
