#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

parser = argparse.ArgumentParser(prog='finalAnnotationAnalysis.py', description='generate lnc/ncRNAs GO frequency and size distribution', add_help=True)
parser.add_argument('-annotation', dest='annotation_file', metavar='guiltByAssociationFunctions.tsv', help='lnc/ncRNAs annotation - guilt by association', type=str, required=True)
parser.add_argument('-obo', dest='obo_file', metavar='go-basic.obo', help='obo file with GO descriptions', type=str, required=True)
parser.add_argument('-out_ncRNA', dest='output_ncRNA', metavar='guiltByAssociation_ncRNAsFunctions.tsv', help='GO ncRNAs functions', type=str, required=True)
parser.add_argument('-out_lncRNA', dest='output_lncRNA', metavar='guiltByAssociation_lncRNAsFunctions.tsv', help='GO lncRNAs functions', type=str, required=True)
parser.add_argument('-plot', dest='plot_file', metavar='guiltByAssociation_lncRNA_size_distribution.png', help='Annotated lncRNAs size distribution', type=str, required=True)
args = parser.parse_args()

annotation_file = args.annotation_file
obo_file = args.obo_file
output_ncRNA = args.output_ncRNA
output_lncRNA = args.output_lncRNA
plot_file = args.plot_file

# Função para contar frequências de termos GO separando lncRNA e ncRNA
def count_go_frequencies(annotation_file):
    ncRNA_go_frequencies = {}
    lncRNA_go_frequencies = {}
    gene_go_seen = set()
    
    df = pd.read_csv(annotation_file, sep='\t')
    
    for index, row in df.iterrows():
        gene = row['Gene']
        go_terms = row['GO (Guilt by association)']
        gene_category = row['Gene Category']
        
        if pd.notna(go_terms):
            go_terms_list = go_terms.split(',')
            
            for term in go_terms_list:
                if (gene, term) not in gene_go_seen:
                    if gene_category == 'lncRNA':
                        if term in lncRNA_go_frequencies:
                            lncRNA_go_frequencies[term] += 1
                        else:
                            lncRNA_go_frequencies[term] = 1
                    elif gene_category == 'ncRNA':
                        if term in ncRNA_go_frequencies:
                            ncRNA_go_frequencies[term] += 1
                        else:
                            ncRNA_go_frequencies[term] = 1
                    gene_go_seen.add((gene, term))  
    
    return lncRNA_go_frequencies, ncRNA_go_frequencies

# Função para ler e processar o arquivo OBO
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

# Função para traduzir os termos GO
def translate_go_terms(go_frequencies, go_terms):
    translated_terms = {}
    for go_id, frequency in go_frequencies.items():
        if go_id in go_terms:
            translated_terms[go_id] = (go_terms[go_id], frequency)
    return translated_terms

# Função para salvar as frequências de termos GO em arquivos
def save_go_frequencies(translated_terms, output_file):
    sorted_terms = sorted(translated_terms.items(), key=lambda item: item[1][1], reverse=True)
    with open(output_file, 'w') as file:
        file.write("GO.ID\tName\tFrequency\n")
        for go_id, (name, frequency) in sorted_terms:
            file.write(f"{go_id}\t{name}\t{frequency}\n")
    print(f"Frequência dos termos GO salva em {output_file}")

# Função para plotar a distribuição do tamanho dos lncRNAs
def plot_lncRNA_transcript_size_distribution(annotation_file):
    df = pd.read_csv(annotation_file, sep='\t')
    
    lncRNA_transcripts = df[df['Transcript Category'] == 'lncRNA']
   
    transcript_sizes = lncRNA_transcripts['Transcript Size']
    transcript_sizes = transcript_sizes.dropna()
    
    plt.figure(figsize=(10, 6))
    sns.histplot(transcript_sizes, color='skyblue')
    
    plt.title('Distribuição do tamanho dos lncRNAs anotados', fontsize=16)
    plt.xlabel('Tamanho do lncRNA', fontsize=14)
    plt.ylabel('Frequência', fontsize=14)
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)

# Contar frequências de GO
lncRNA_go_frequencies, ncRNA_go_frequencies = count_go_frequencies(annotation_file)
go_terms = parse_go_obo(obo_file)

# Traduzir e salvar frequências de termos GO para os lncRNAs
translated_lncRNA_terms = translate_go_terms(lncRNA_go_frequencies, go_terms)
save_go_frequencies(translated_lncRNA_terms, output_lncRNA)

# Traduzir e salvar frequências de termos GO para os ncRNAs
translated_ncRNA_terms = translate_go_terms(ncRNA_go_frequencies, go_terms)
save_go_frequencies(translated_ncRNA_terms, output_ncRNA)

# Plotr distribuição do tamanho dos lncRNAs
plot_lncRNA_transcript_size_distribution(annotation_file)
