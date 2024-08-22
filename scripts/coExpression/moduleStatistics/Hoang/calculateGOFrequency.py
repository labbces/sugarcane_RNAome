#!/usr/bin/env python

import pandas as pd
import os
import csv

# Caminho para os arquivos
module_file = '../moduleCNCwithlncRNAs.tsv'
go_directory = './'
obo_file = 'go-basic.obo'
output_file = 'GOFrequencies.tsv'

def calculate_go_frequencies(module_file, go_directory):
    modules_df = pd.read_csv(module_file, sep='\t')
    go_frequencies = {}

    for module in modules_df['Module']:
        go_file = os.path.join(go_directory, f'module_{module}.csv')
        
        if os.path.exists(go_file):
            with open(go_file, 'r') as file:
                header = file.readline().strip().split(',')
                go_id_index = header.index('GO.ID')
                    
                for line in file:
                    parts = line.strip().split(',')
                    if len(parts) > go_id_index:
                        go_id = parts[go_id_index].strip()
                        if go_id:
                            if go_id in go_frequencies:
                                go_frequencies[go_id] += 1
                            else:
                                go_frequencies[go_id] = 1

    
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
