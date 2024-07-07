#!/usr/bin/env python

import pandas as pd
import numpy as np

input_file = 'updated_panRNAome_Length.tsv'
df = pd.read_csv(input_file, sep='\t', header=None)

df.columns = ['panRNAome_category', 'gene_name', 'transcript_name', 'gene_function', 'transcript_length', 'transcript_function']

# Ordem de precedência das funções
function_precedence = {
    'protein and non-coding': 1,
    #'protein and lncRNA': 2,
    #'protein-coding': 3,
    'ncRNA': 2,
    'lncRNA': 3
}

# Função para determinar a função do gene com base nos transcritos
def determine_gene_function(transcript_functions):
    # Filtra valores não nulos
    valid_functions = set(transcript_functions.dropna())
    
    # Verifica se há transcritos "protein-coding" e "ncRNA"
    if 'protein-coding' in valid_functions and 'ncRNA' in valid_functions:
        return 'protein and non-coding'
    
    # Ordena as funções de acordo com a precedência
    sorted_functions = sorted(valid_functions, key=lambda x: function_precedence.get(x, float('inf')))
    return sorted_functions[0] if sorted_functions else np.nan  # Retorna a função de maior precedência ou NA se estiver vazio

# Agrupamento por gene e determinação da função do gene
gene_functions = df.groupby('gene_name')['transcript_function'].apply(determine_gene_function).reset_index()
gene_functions.columns = ['gene_name', 'updated_gene_function']

# Atualização da coluna de função do gene no dataframe original
df = df.merge(gene_functions, on='gene_name', how='left')
df['gene_function'] = df['updated_gene_function']
df.drop(columns=['updated_gene_function'], inplace=True)

# Convertendo transcript_length para inteiro (removendo ".0" se for inteiro)
df['transcript_length'] = df['transcript_length'].astype('Int64')  # Utiliza tipo Int64 para manter NA e inteiros

output_file = 'updated_panRNAome_GeneFunction_Length.tsv'
df.to_csv(output_file, sep='\t', header=False, index=False, na_rep='NA')
