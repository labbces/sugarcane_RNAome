#!/usr/bin/env python

import pandas as pd
import argparse
import sys

parser = argparse.ArgumentParser(description='compute module stats')
parser.add_argument('-m', metavar='cliques.csv', dest='cliques', type=str, help='file with cliques (modules)', required=True)
parser.add_argument('-f', metavar='filtered_cliques.csv', dest='filtered_cliques', type=str, help='file with filtered cliques (cliques with more than 5 genes)', required=True)
parser.add_argument('-a', metavar='annotation.tsv', dest='annotation', type=str, help='file with genes annotation', required=True)
parser.add_argument('-o', metavar='output.tsv', dest='output', type=str, help='output file with modules statistics', required=True)
args = parser.parse_args()

cliques = args.cliques
filtered_cliques = args.filtered_cliques
annotation = args.annotation
output = args.output

A = pd.read_csv(cliques, header=None, names=["Gene", "Pertinency", "Module"], sep=' ')
B = pd.read_csv(filtered_cliques, header=None, names=["Module"], sep=' ')

dtype_dict = {
    "Category": str,
    "Gene": str,
    "Transcript": str,
    "Gene Category": str,
    "Transcript Size": str,
    "Transcript Category": str,
    "Transcript Rfam family": str
}

C = pd.read_csv(annotation, sep="\t", dtype=dtype_dict, index_col=False)#, names=["Category", "Gene", "Transcript", "Gene Category", "Transcript Size", "Transcript Category", "Transcript Rfam family"], dtype=dtype_dict)

#print(A)
#print(B)
print(C.head())
print(C.columns)
print(C.dtypes)

#sys.exit()

# filtrar modulos de interesse
A_filtered = A[A['Module'].isin(B['Module'])]
#print(A_filtered)

# agrupar por modulo e contar as categorias de genes (protein, protein and non-coding, ncRNA, lncRNA)
result = A_filtered.merge(C[['Gene', 'Gene Category']], on='Gene')
summary = result.groupby('Module').agg(
    Genes=('Gene', 'count'),
    Genes_coding=('Gene Category', lambda x: (x == 'protein-coding').sum()),
    Genes_protein_and_non_coding=('Gene Category', lambda x: (x == 'protein and non-coding').sum()),
    Genes_ncRNA=('Gene Category', lambda x: (x == 'ncRNA').sum()),
    Genes_lncRNA=('Gene Category', lambda x: (x == 'lncRNA').sum())
).reset_index()
#print(summary)

summary.to_csv(output, sep='\t',index=False)
