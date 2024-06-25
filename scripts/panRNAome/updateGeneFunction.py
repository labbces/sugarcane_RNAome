#!/usr/bin/env python

import pandas as pd

input_file = 'updated_panTranscriptome_panRNAome_Length.tsv'
df = pd.read_csv(input_file, sep='\t', header=None)

df.columns = ['panRNAome_category', 'gene_name', 'transcript_name', 'gene_function', 'transcript_length', 'transcript_function']

function_precedence = {
    'protein and non-coding': 1,
    'protein and lncRNA': 2,
    'protein-coding': 3,
    'ncRNA': 4,
    'lncRNA': 5,
}

def determine_gene_function(transcript_functions):
    # Ordena as funções de acordo com a precedência
    sorted_functions = sorted(transcript_functions, key=lambda x: function_precedence.get(x, float('inf')))
    return sorted_functions[0]  # Retorna a função de maior precedência

# Agrupamento por gene e determinação da função do gene
gene_functions = df.groupby('gene_name')['transcript_function'].apply(determine_gene_function).reset_index()
gene_functions.columns = ['gene_name', 'updated_gene_function']

df = df.drop(columns=['gene_function']).merge(gene_functions, on='gene_name')
df = df.rename(columns={'updated_gene_function': 'gene_function'})

df = df[['panRNAome_category', 'gene_name', 'transcript_name', 'gene_function', 'transcript_length', 'transcript_function']]

output_file = 'updated_panTranscriptome_panRNAome_GeneFunction_Length.tsv'
df.to_csv(output_file, sep='\t', header=False, index=False)
