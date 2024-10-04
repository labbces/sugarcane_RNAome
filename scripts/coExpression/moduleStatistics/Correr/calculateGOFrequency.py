#!/usr/bin/env python

import pandas as pd
import os
import csv
import argparse

parser = argparse.ArgumentParser(prog='calculateGOFrequency.py', description='Calculate GO Frequency of enriched CNC modules with lncRNAs', add_help=True)
parser.add_argument('-modules', dest='module_file', metavar='moduleCNCwithlncRNAs.tsv', help='Table containing modules CNC with lncRNAs', type=str, required=True)
parser.add_argument('-GO', dest='go_directory', metavar='./', help='Directory containing GO terms of enriched modules', type=str, required=True)
parser.add_argument('-obo', dest='obo_file', metavar='go-basic.obo', help='Gene Ontology descriptions', type=str, required=True)
parser.add_argument('-out_freq', dest='output_file', metavar='GOFrequencies.tsv', help='Output file with GO Frequency', type=str, required=True)
parser.add_argument('-out_modules', dest='output_modulos_enriquecidos', metavar='moduleCNCwithlncRNAswithEnrichedGO.tsv', help='Table with enriched modules', type=str, required=True)
args = parser.parse_args()

# Caminho para os arquivos
#module_file = '../moduleCNCwithlncRNAs.tsv'
#go_directory = './'
#obo_file = 'go-basic.obo'
#output_file = 'GOFrequencies.tsv'
#output_modulos_enriquecidos = 'moduleCNCwithlncRNAswithEnrichedGO.tsv'

module_file = args.module_file
go_directory = args.go_directory
obo_file = args.obo_file
output_file = args.output_file
output_modulos_enriquecidos = args.output_modulos_enriquecidos


def calculate_go_frequencies(module_file, go_directory):
    modules_df = pd.read_csv(module_file, sep='\t')
    go_frequencies = {}
    
    total_genes = 0
    total_genes_coding = 0
    total_genes_ncRNA = 0
    total_genes_lncRNA = 0
    
    modules_with_enriched_go_terms = []
    
    for index, row in modules_df.iterrows():
        module = row['Module']
        go_file = os.path.join(go_directory, f'module_{module}.csv')
        
        if os.path.exists(go_file):
            try:
                go_df = pd.read_csv(go_file, sep=',', usecols=[0])
                
                if not go_df.empty and 'GO.ID' in go_df.columns and len(go_df) > 0:
                    total_genes += row['Genes']
                    total_genes_coding += row['Genes_coding']
                    total_genes_ncRNA += row['Genes_ncRNA']
                    total_genes_lncRNA += row['Genes_lncRNA']
                    
                    modules_with_enriched_go_terms.append(module)
                    
                    for go_id in go_df['GO.ID'].dropna():
                        go_id = go_id.strip()
                        if go_id:
                            if go_id in go_frequencies:
                                go_frequencies[go_id] += 1
                            else:
                                go_frequencies[go_id] = 1
            except pd.errors.ParserError as e:
                print(f"Erro ao ler o arquivo {go_file} (Módulo {module}): {e}")
    
    print(f"Modulos com termos enriquecidos: {modules_with_enriched_go_terms}")
    print(f"Total de Genes: {total_genes}")
    print(f"Total de Genes Coding: {total_genes_coding}")
    print(f"Total de Genes ncRNA: {total_genes_ncRNA}")
    print(f"Total de Genes lncRNA: {total_genes_lncRNA}")
    
    enriched_modules_df = pd.DataFrame(modules_with_enriched_go_terms, columns=['Modulos com termos enriquecidos'])
    enriched_modules_df.to_csv(output_modulos_enriquecidos, index=False, sep='\t')
        
    return go_frequencies

def parse_go_obo(obo_file):
    go_terms = {}
    with open(obo_file, 'r') as file:
        term_id = None
        for line in file:
            if line.startswith('[Term]'):
                term_id = None
            elif line.startswith('id:'):
                term_id = line.strip().split(': ')[1]
            elif line.startswith('name:') and term_id:
                go_terms[term_id] = line.strip().split(': ')[1]
                term_id = None
    return go_terms

def translate_go_terms(go_frequencies, go_terms):
    translated_terms = {}
    for go_id, frequency in go_frequencies.items():
        if go_id in go_terms:
            translated_terms[go_id] = (go_terms[go_id], frequency)
    return translated_terms

def save_frequencies(translated_terms, output_file):
    df = pd.DataFrame.from_dict(translated_terms, orient='index', columns=['Name', 'Frequency'])
    df.index.name = 'GO.ID'
    df = df.reset_index()
    df = df.sort_values(by='Frequency', ascending=False)
    df.to_csv(output_file, index=False, quoting=csv.QUOTE_NONE, sep='\t')

go_frequencies = calculate_go_frequencies(module_file, go_directory)

go_terms = parse_go_obo(obo_file)

translated_terms = translate_go_terms(go_frequencies, go_terms)

save_frequencies(translated_terms, output_file)
