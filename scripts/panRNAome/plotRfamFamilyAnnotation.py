#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

annotation_file = 'updated_panRNAome_GeneFunction_Length_with_Rfam.tsv'
columns = ['panRNAome_category', 'gene_name', 'transcript_name', 'gene_function', 'transcript_length', 'transcript_function', 'Transcript Rfam family', 'Gene Rfam family']
annotation_df = pd.read_csv(annotation_file, sep='\t', header=None, names=columns)

# nova coluna com presença ou ausencia de Gene Rfam family
annotation_df['Has Transcript Rfam family'] = annotation_df['Transcript Rfam family'].notna()

# agrupar por gene_name e gene_function e determinar se o gene possui uma "Gene Rfam family"
#gene_df = annotation_df.groupby(['gene_name', 'gene_function']).agg({'Has Gene Rfam family': 'any'}).reset_index()
gene_df = annotation_df.groupby(['transcript_name', 'transcript_function']).agg({'Has Transcript Rfam family': 'any'}).reset_index()

# plot quantidade de genes em cada classificação usando a coluna 'gene_function'
plt.figure(figsize=(10, 6))
ax = sns.countplot(data=gene_df, x='transcript_function')
plt.title('Distribuição das categorias dos transcritos')
plt.xlabel('Categoria dos transcritos')
plt.ylabel('Frequência dos transcritos')
plt.xticks(rotation=45)
plt.tight_layout()

# add label de contagem acima de cada barra
for p in ax.patches:
    ax.annotate(format(p.get_height(), '.0f'), 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha='center', va='bottom')

plt.savefig('transcript_function_distribution.png')
#plt.show()
plt.clf()  # Limpar a figura atual

# filtrar genes com gene_function (lncRNA, protein and lncRNA, ncRNA, protein and non-coding)
#filtered_gene_df = gene_df.query('gene_function == "lncRNA" or gene_function == "protein and lncRNA" or gene_function == "ncRNA" or gene_function == "protein and non-coding"')
filtered_gene_df = gene_df.query('transcript_function == "lncRNA" or transcript_function == "ncRNA" or transcript_function == "protein and non-coding"')

# plot quantidade de genes com e sem "Gene Rfam family"
plt.figure(figsize=(14, 8))

##ax = sns.countplot(data=filtered_gene_df, x='transcript_function', hue='Has Transcript Rfam family')
ax = sns.countplot(data=gene_df, x='transcript_function', hue='Has Transcript Rfam family')

plt.title('Distribuição de ncRNAs classificados pelo Rfam')
plt.xlabel('Categoria do grupo')
plt.ylabel('Frequência de transcritos')
plt.legend(title='Transcrito classificado pelo Rfam', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=45)
plt.tight_layout()

# add label de contagem acima de cada barra
for p in ax.patches:
    ax.annotate(format(p.get_height(), '.0f'),
                (p.get_x() + p.get_width() / 2., p.get_height()),
                ha='center', va='bottom')

plt.savefig('transcript_function_by_Rfam_family_presence.png')
#plt.show()
plt.clf()  # Limpar a figura atual
