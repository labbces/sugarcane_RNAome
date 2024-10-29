#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt

df = 'updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv'
columns = ['Category', 'Gene', 'Transcript', 'Gene Category', 'Transcript Size', 'Transcript Category', 'Transcript Rfam family', 'Transcript GO', 'Gene GO']
annotation_df = pd.read_csv(df, sep='\t', header=None, names=columns)

df_exclusive = annotation_df[annotation_df['Category'] == 'Exclusive']

df_exclusive['Genotype'] = df_exclusive['Transcript'].str.split('_').str[0]

df_unique_genes = df_exclusive[['Genotype', 'Gene', 'Gene Category']].drop_duplicates()

gene_counts = df_unique_genes.groupby(['Genotype', 'Gene Category']).size().unstack(fill_value=0)

plt.figure(figsize=(14, 8))
colors = plt.get_cmap("tab20").colors
gene_counts.plot(kind='bar', stacked=False, color=colors, edgecolor='black', ax=plt.gca())

plt.xlabel('Genótipo')
plt.ylabel('Número de genes exclusivos')
plt.xticks(rotation=90)
plt.legend(title='Categoria do gene', bbox_to_anchor=(1.05, 1), loc='upper left') 
plt.tight_layout()

plt.savefig('exclusiveGenesPerGenotype_by_GeneCategory.png')
#plt.show()
plt.clf()
