#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(prog='plotCategoryOfGenesPerGenotype.py', description='plot pan category of genes per genotype', add_help=True)
parser.add_argument('-a', '--annotation', dest='df', metavar='updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv',help="panTranscriptome and panRNAome annotation", required=True)
parser.add_argument('-c', '--category', dest='category', metavar='Exclusive', help="pan category", required=True)

args = parser.parse_args()
df = args.df
category = args.category
output_category = category.lower()

columns = ['Category', 'Gene', 'Transcript', 'Gene Category', 'Transcript Size', 'Transcript Category', 'Transcript Rfam family', 'Transcript GO', 'Gene GO']
annotation_df = pd.read_csv(df, sep='\t', header=None, names=columns)

df_category = annotation_df[annotation_df['Category'] == category]

df_category['Genotype'] = df_category['Transcript'].str.split('_').str[0]

df_unique_genes = df_category[['Genotype', 'Gene', 'Gene Category']].drop_duplicates()

gene_counts = df_unique_genes.groupby(['Genotype', 'Gene Category']).size().unstack(fill_value=0)

plt.figure(figsize=(14, 8))
colors = plt.get_cmap("tab20").colors
gene_counts.plot(kind='bar', stacked=False, color=colors, edgecolor='black', ax=plt.gca())

plt.xlabel('Genótipo')
plt.ylabel('Número de genes' + output_category)
plt.xticks(rotation=90)
plt.legend(title='Categoria do gene', bbox_to_anchor=(1.05, 1), loc='upper left') 
plt.tight_layout()

plt.savefig(output_category + 'GenesPerGenotype_by_GeneCategory.png')
#plt.show()
plt.clf()
