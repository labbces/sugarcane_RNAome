#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(prog='getGOTermsUsedToAnnotateEnrichedModules.py', description='Get GO terms from proteins used to annotate CNC modules with lncRNAs', add_help=True)
parser.add_argument('-modules', dest='modules_file', metavar='moduleCNCwithlncRNAswithEnrichedGO.tsv', help='Table with enriched modules', type=str, required=True)
parser.add_argument('-genes', dest='genes_modules_file', metavar='cliques.csv', help='Table with genes per module', type=str, required=True)
parser.add_argument('-annotation', dest='gene_annotation_file', metavar='updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv', help='Annotation table', type=str, required=True)
parser.add_argument('-out', dest='output_file', metavar='GOTermsUsedInModuleCNCwithlncRNAswithEnrichedGO.tsv', help='Output table with GO per genes in enriched CNC modules', type=str, required=True)
args = parser.parse_args()

#modules_file = '../enrichedGOFrequency/moduleCNCwithlncRNAswithEnrichedGO.tsv'
#genes_modules_file = 'Correr2020_counts_filters_VST_topCV_mcl_formated_cliques.csv'
#gene_annotation_file = 'updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv'
#output_file = 'GOTermsUsedInModuleCNCwithlncRNAswithEnrichedGO.tsv'

modules_file = args.modules_file
genes_modules_file = args.genes_modules_file
gene_annotation_file = args.gene_annotation_file
output_file = args.output_file

def count_go_terms(modules_file, genes_modules_file, gene_annotation_file, output_file):
    modules_with_enriched_go = pd.read_csv(modules_file, sep='\t')

    genes_modules_df = pd.read_csv(genes_modules_file, sep=' ', header=None, names=['Gene', 'Pertinency', 'Module'])

    #gene_annotation_df = pd.read_csv(gene_annotation_file, sep='\t', index_col=False, usecols=['Gene', 'Gene GO'])
    gene_annotation_df = pd.read_csv(gene_annotation_file, sep='\t', index_col=False, usecols=['Gene', 'Transcript GO'])

    go_terms_by_module = {}

    for module in modules_with_enriched_go['Modulos com termos enriquecidos']:

        genes_in_module = genes_modules_df[genes_modules_df['Module'] == module]['Gene'].tolist()

        #go_terms = gene_annotation_df[gene_annotation_df['Gene'].isin(genes_in_module)]['Gene GO'].dropna().unique()
        go_terms = gene_annotation_df[gene_annotation_df['Gene'].isin(genes_in_module)]['Transcript GO'].dropna().unique()
        
        go_terms_by_module[module] = ','.join(go_terms)
        
    go_terms_df = pd.DataFrame({
        'Modulos com termos enriquecidos': go_terms_by_module.keys(),
        'Termos GO utilizados para enriquecimento neste modulo': go_terms_by_module.values()
    })

    go_terms_df.to_csv(output_file, sep='\t', index=False)

    print(f"Tabela de termos GO por módulo salva em {output_file}")

count_go_terms(modules_file, genes_modules_file, gene_annotation_file, output_file)
