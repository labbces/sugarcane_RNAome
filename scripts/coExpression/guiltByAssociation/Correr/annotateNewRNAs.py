#!/usr/bin/env python

import argparse
import os
import csv
import pandas as pd

parser = argparse.ArgumentParser(description='transfer GO to ncRNAs - Guilt by association')
parser.add_argument('-m', metavar='cliques.csv', dest='cliques', type=str, help='file with cliques (modules)', required=True)
parser.add_argument('-go', metavar='enrichedGOFrequency', dest='enrichedGO', type=str, help='directory with enriched GO terms for each module', required=True)
parser.add_argument('-a', metavar='annotation.tsv', dest='annotation', type=str, help='file with genes annotation', required=True)
parser.add_argument('-o', metavar='output.tsv', dest='output', type=str, help='output file with modules statistics', required=True)
args = parser.parse_args()

cliques = args.cliques
enrichedGO = args.enrichedGO
annotation = args.annotation
output = args.output

def extrair_termos_GO(directory):
    go_terms = {}
    for filename in os.listdir(directory):
        if filename.startswith('module_'):
            module_num = filename.split('_')[1].split('.')[0] 
            with open(os.path.join(directory, filename)) as csv_file:
                reader = csv.DictReader(csv_file)
                for row in reader:
                    if module_num not in go_terms:
                        go_terms[module_num] = []
                    go_terms[module_num].append(row['GO.ID'])
    #print(go_terms)
    return go_terms
	
def obter_genes_modulos(filepath):
    genes_por_modulo = {}
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split()
            gene, modulo = parts[0], parts[2]
            if modulo not in genes_por_modulo:
                genes_por_modulo[modulo] = []
            genes_por_modulo[modulo].append(gene)
    #print(genes_por_modulo)
    return genes_por_modulo

def filtrar_genes_ncRNAs(filepath):
    df = pd.read_csv(filepath, sep='\t', index_col=False, usecols=['Gene', 'Gene Category']) 
    #print(df)
    df['Gene Category'] = df['Gene Category'].astype(str).fillna('').str.strip()

    ncRNA_genes = df[df['Gene Category'].str.contains('ncRNA|lncRNA', na=False)]

    ncRNA_genes = ncRNA_genes.drop_duplicates(subset='Gene')

    #print(ncRNA_genes)
    
    return ncRNA_genes

def atribuir_go_terms(gene, genes_por_modulo, go_terms):
    for modulo, genes in genes_por_modulo.items():
        if gene in genes:
            return ','.join(go_terms.get(modulo, []))
    return ''

def adicionar_termos_GO(ncRNA_df, genes_por_modulo, go_terms):
    ncRNA_df['GO (Guilt by association)'] = ncRNA_df['Gene'].apply(lambda gene: atribuir_go_terms(gene, genes_por_modulo, go_terms))
    # print(ncRNA_df)
    return ncRNA_df

def gerar_tabela_final(enriched_go_dir, genes_modulos_filepath, annot_filepath, output_filepath):
    termos_go = extrair_termos_GO(enriched_go_dir)
    
    genes_por_modulo = obter_genes_modulos(genes_modulos_filepath)
    
    ncRNAs = filtrar_genes_ncRNAs(annot_filepath)
    
    tabela_final = adicionar_termos_GO(ncRNAs, genes_por_modulo, termos_go)
    
    tabela_final.to_csv(output_filepath, sep='\t', index=False)

gerar_tabela_final(enrichedGO, cliques, annotation, output)
