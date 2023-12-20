#!/usr/bin/env python3

def load_table(file_path):
    data = set()
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            if len(columns) >= 3:
                data.add(columns[2])
    return data

def find_exact_matches(table1, table2):
    exact_matches = table1.intersection(table2)
    return exact_matches

file_path1 = "panTranscriptomeClassificationTable_I6.0.tsv"
file_path2 = "panRNAomeClassificationTable_0.8.tsv"

# load tables
table1 = load_table(file_path1)
table2 = load_table(file_path2)

exact_matches = find_exact_matches(table1, table2)
  
output_file_path = "matches_panTranscriptome_panRNAome.tsv"

with open(output_file_path, 'w') as output_file:
    output_file.write("sequence matches in panTranscriptomeClassificationTable_I6.0.tsv and panRNAomeClassificationTable_0.8.tsv:\n")
    for match in exact_matches:
        output_file.write(match + '\n')

print(f"Exact matches saved in: {output_file_path}")
