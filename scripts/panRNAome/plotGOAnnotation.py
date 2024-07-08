#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

annotation_file = 'updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv'
columns = ['panRNAome_category', 'gene_name', 'transcript_name', 'gene_function', 'transcript_length', 'transcript_function', 'Transcript Rfam family', 'Transcript GO','Gene GO']
annotation_df = pd.read_csv(annotation_file, sep='\t', header=None, names=columns)

# nova coluna com presença ou ausencia de Gene Rfam family
annotation_df['Has Gene GO'] = annotation_df['Gene GO'].notna()
annotation_df['Has Transcript GO'] = annotation_df['Transcript GO'].notna()

# agrupar por gene_name e gene_function e determinar se o gene possui uma "Gene Rfam family"
gene_df = annotation_df.groupby(['gene_name', 'gene_function']).agg({'Has Gene GO': 'any'}).reset_index()
transcript_df = annotation_df.groupby(['transcript_name', 'transcript_function']).agg({'Has Transcript GO': 'any'}).reset_index()

# plot quantidade de genes em cada classificação usando a coluna 'gene_function'
plt.figure(figsize=(10, 6))
ax = sns.countplot(data=gene_df, x='gene_function')
plt.title('Distribuição de grupos com Gene Ontology (BP)')
plt.xlabel('Categoria dos grupos')
plt.ylabel('Frequência dos grupos')
plt.xticks(rotation=45)
plt.tight_layout()

# add label de contagem acima de cada barra
for p in ax.patches:
    ax.annotate(format(p.get_height(), '.0f'), 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha='center', va='bottom')

plt.savefig('gene_GO_distribution.png')
#plt.show()
plt.clf()  # Limpar a figura atual

# filtrar genes com gene_function (lncRNA, protein and lncRNA, ncRNA, protein and non-coding)
#filtered_gene_df = gene_df.query('gene_function == "lncRNA" or gene_function == "ncRNA" or gene_function == "protein-coding" or gene_function == "protein and non-coding"')

# plot quantidade de genes com e sem "Gene GO"
plt.figure(figsize=(14, 8))
ax = sns.countplot(data=gene_df, x='gene_function', hue='Has Gene GO')
plt.title('Distribuição de grupos com Gene Ontology (BP)')
plt.xlabel('Categoria dos grupos')
plt.ylabel('Frequência de grupos')
plt.legend(title='Grupo classificado com termo GO', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=45)
plt.tight_layout()

# add label de contagem acima de cada barra
for p in ax.patches:
    ax.annotate(format(p.get_height(), '.0f'),
                (p.get_x() + p.get_width() / 2., p.get_height()),
                ha='center', va='bottom')

plt.savefig('gene_function_by_GO_presence.png')
#plt.show()
plt.clf()  # Limpar a figura atual

# plot quantidade de transcritos com e sem "Transcript GO"
plt.figure(figsize=(14, 8))
ax = sns.countplot(data=transcript_df, x='transcript_function', hue='Has Transcript GO')
plt.title('Distribuição de transcritos com Gene Ontology (BP)')
plt.xlabel('Categoria dos transcritos')
plt.ylabel('Frequência de transcritos')
plt.legend(title='Transcrito classificado com termo GO', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=45)
plt.tight_layout()

# add label de contagem acima de cada barra
for p in ax.patches:
    ax.annotate(format(p.get_height(), '.0f'),
                (p.get_x() + p.get_width() / 2., p.get_height()),
                ha='center', va='bottom')

plt.savefig('transcript_function_by_GO_presence.png')
#plt.show()
plt.clf()  # Limpar a figura atual
