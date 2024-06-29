#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

file_path = 'updated_panTranscriptome_panRNAome_GeneFunction_Length.tsv'
data = pd.read_csv(file_path, sep='\t', header=None, names=['panRNAome category', 'Gene', 'Transcript', 'Gene Function', 'Transcript Size', 'Transcript Function'])

# Extrair o genótipo do nome do transcrito
data['Genotype'] = data['Transcript'].apply(lambda x: x.split('_')[0])

# Remover linhas com valores NA em 'Transcript Size' e 'Transcript Function'
data = data.dropna(subset=['Transcript Size', 'Transcript Function'])
data['Transcript Size'] = data['Transcript Size'].astype(int)

# Transformar Transcript Size em inteiro (converter NAs para 0)
#data['Transcript Size'] = data['Transcript Size'].fillna(0).astype(int)

# 1) Quantos transcritos tem os genes? Qual a distribuição do tamanho deles?
transcripts_per_gene = data.groupby('Gene')['Transcript'].count()
transcript_size_distribution = data['Transcript Size']
transcripts_per_gene_function = data.groupby('Gene Function')['Transcript'].count()

print(f'Número de transcritos por gene: {transcripts_per_gene}')
print(f'Distribuição do tamanho dos transcritos: {transcript_size_distribution}')
print(f'Número de transcritos por funçao do gene: {transcripts_per_gene_function}')

# Histograma de número de transcritos por gene
plt.figure(figsize=(14, 6))
plt.subplot(1, 2, 1)
sns.histplot(transcripts_per_gene, bins=range(0, transcripts_per_gene.max() + 10, 10), kde=True)
plt.title('Distribuição do número de transcritos por gene')
plt.xlabel('Número de transcritos')
plt.ylabel('Frequência')

# Histograma de número de transcritos por função do gene
plt.subplot(1, 2, 2)
colors = sns.color_palette('husl', len(transcripts_per_gene_function))

for i, (gene_function, count) in enumerate(transcripts_per_gene_function.items()):
    sns.histplot([count], bins=range(0, transcripts_per_gene_function.max() + 10, 10), color=colors[i], label=gene_function, kde=True)

plt.title('Distribuição do número de transcritos por função do gene')
plt.xlabel('Número de transcritos')
plt.ylabel('Frequência')
plt.legend(title='Função do gene', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
# plt.show()
plt.savefig('transcripts_distribution.png')
plt.clf()  # Limpar a figura atual

plt.figure(figsize=(10, 6))
sns.histplot(transcript_size_distribution, kde=True, bins=range(0, transcript_size_distribution.max() + 10, 10))
plt.title('Distribuição do tamanho dos transcritos')
plt.xlabel('Tamanho do transcrito')
plt.ylabel('Frequência')
plt.tight_layout()
# plt.show()
plt.savefig('transcript_size_distribution.png')
plt.clf()  # Limpar a figura atual

# 2) Quantos transcritos são lncRNAs?
lncRNA_transcripts = data[data['Transcript Function'].str.contains('lncRNA')]
num_lncRNA_transcripts = len(lncRNA_transcripts)

print(f'Número de transcritos lncRNAs: {num_lncRNA_transcripts}')

# 3) Nos genes com transcritos ncRNAs, como estão os transcritos?
ncRNA_categories = ['ncRNA', 'lncRNA', 'protein and non-coding', 'protein and lncRNA']
ncRNA_transcripts = data[data['Transcript Function'].isin(ncRNA_categories)]

plt.figure(figsize=(14, 6))

# Histograma de tamanho dos transcritos ncRNAs
plt.subplot(1, 2, 1)
sns.histplot(ncRNA_transcripts['Transcript Size'], kde=True, bins=range(0, ncRNA_transcripts['Transcript Size'].max() + 10, 10))
plt.title('Distribuição do tamanho dos ncRNAs')
plt.xlabel('Tamanho do transcrito')
plt.ylabel('Frequência')

# Histograma colorido por categorias de ncRNA
plt.subplot(1, 2, 2)
colors = sns.color_palette('husl', len(ncRNA_categories))

for i, category in enumerate(ncRNA_categories):
    category_data = ncRNA_transcripts[ncRNA_transcripts['Transcript Function'] == category]
    sns.histplot(category_data['Transcript Size'], kde=True, color=colors[i], label=category, bins=range(0, ncRNA_transcripts['Transcript Size'].max() + 10, 10))

plt.title('Distribuição do tamanho dos ncRNAs')
plt.xlabel('Tamanho do transcrito')
plt.ylabel('Frequência')
plt.legend(title='Categoria de ncRNA', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
# plt.show()
plt.savefig('ncRNA_transcripts_size_distribution.png')
plt.clf()  # Limpar a figura atual

# 4) Qual a distribuição dos genótipos nos genes lncRNAs? Todos os 48 genótipos estão presentes?
lncRNA_genes = lncRNA_transcripts['Gene'].unique()
genotypes_in_lncRNA_genes = data[data['Gene'].isin(lncRNA_genes)]['Genotype'].unique()

print(f'Número de genótipos nos genes lncRNAs: {len(genotypes_in_lncRNA_genes)}')
print(f'Genótipos presentes nos genes lncRNAs: {genotypes_in_lncRNA_genes}')

# 5) Distribuição dos tamanhos dos transcritos por genótipo
plt.figure(figsize=(20, 6))  
sns.boxplot(x='Genotype', y='Transcript Size', data=data)
plt.xticks(rotation=90)
plt.title('Distribuição do tamanho dos transcritos por genótipo')
plt.xlabel('Genótipo')
plt.ylabel('Tamanho do transcrito')
plt.tight_layout()
# plt.show()
plt.savefig('transcript_size_by_genotype.png')
plt.clf()  # Limpar a figura atual

# 6) Comparação das funções dos transcritos entre os genótipos
plt.figure(figsize=(20, 8))  
sns.countplot(x='Genotype', hue='Transcript Function', data=data)
plt.xticks(rotation=90)
plt.title('Funções dos transcritos por genótipo')
plt.xlabel('Genótipo')
plt.ylabel('Frequência')
plt.tight_layout()
# plt.show()
plt.savefig('transcript_function_by_genotype.png')
plt.clf()  # Limpar a figura atual

# 7) Nos genes classificados como “ncRNA”, “lncRNA”, “protein and non-coding” e “protein and lncRNA” como é a distribuição dos transcritos? 
# Tem transcritos ncRNA junto com lncRNAs? São apenas lncRNAs juntos? Tem transcritos “protein” e “lncRNA” juntos?

genes_interesse = data[data['Gene Function'].isin(ncRNA_categories)]
transcript_function_by_gene_function = genes_interesse.groupby(['Gene Function', 'Transcript Function']).size().unstack(fill_value=0)

fig, ax = plt.subplots(figsize=(10, 6))
transcript_function_by_gene_function.plot(kind='bar', stacked=True, ax=ax)
ax.set_title('Distribuição dos transcritos por função dos genes')
ax.set_xlabel('Gene function')
ax.set_ylabel('Frequência de transcritos')
plt.xticks(rotation=45)
plt.tight_layout()
plt.legend(title='Transcript function')
#plt.show()
plt.savefig('transcript_function_by_gene_function.png')
plt.clf()  # Limpar a figura atual
