#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

file_path = 'updated_panTranscriptome_GeneFunction_Length_with_Rfam_with_GO.tsv'

# Ler o CSV sem definir uma coluna como índice
data = pd.read_csv(file_path, sep='\t', header=None, names=[
    'panRNAome category', 'Gene', 'Transcript', 'Gene Function', 
    'Transcript Size', 'Transcript Function', 'Transcript Rfam family', 
    'Gene Rfam family'], index_col=False, low_memory=False)

# Transformar Transcript em str (apenas para panTranscriptoma)
data['Transcript'] = data['Transcript'].astype(str)

# Remover valores ausentes na coluna 'Transcript Function'
data = data.dropna(subset=['Transcript Function'])

# Redefinir o índice para garantir que ele seja único
data = data.reset_index(drop=True)

# Extrair o genótipo do nome do transcrito
data['Genotype'] = data['Transcript'].apply(lambda x: x.split('_')[0])

# Remover linhas com valores NA em 'Transcript Size'
data = data.dropna(subset=['Transcript Size'])

# Função para verificar se um valor é numérico
def is_numeric(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# Transformar Transcript Size em inteiro, substituindo valores não numéricos por 0
data['Transcript Size'] = data['Transcript Size'].apply(lambda x: int(float(x)) if is_numeric(x) else 0)

# Verificar se a coluna 'Transcript Size' não está vazia
if data['Transcript Size'].empty:
    raise ValueError("A coluna 'Transcript Size' está vazia após a limpeza dos dados.")





# 1) Quantos transcritos tem os genes? Qual a distribuição do tamanho deles?
transcripts_per_gene = data.groupby('Gene')['Transcript'].count()
transcript_size_distribution = data['Transcript Size']
transcripts_per_gene_function = data.groupby('Gene Function')['Transcript'].count()

print(f'Número de transcritos por grupo: {transcripts_per_gene}')
print(f'Distribuição do tamanho dos transcritos: {transcript_size_distribution}')
print(f'Número de transcritos por categoria do grupo: {transcripts_per_gene_function}')

# Histograma de número de transcritos por gene

# Calcular a média e a mediana
mean_transcripts = transcripts_per_gene.mean()
median_transcripts = transcripts_per_gene.median()

plt.figure(figsize=(10, 6))
sns.histplot(transcripts_per_gene, bins=range(0, transcripts_per_gene.max() + 10, 10), kde=False)
plt.title('Distribuição do número de transcritos por grupo')
plt.xlabel('Número de transcritos')
plt.ylabel('Frequência')
plt.yscale('log') 
plt.xticks(rotation=45, fontsize = 'x-small')

# Adicionar linhas verticais para média e mediana
plt.axvline(mean_transcripts, color='r', linestyle='--', linewidth=2, label=f'Média: {mean_transcripts:.2f}')
plt.axvline(median_transcripts, color='b', linestyle='-', linewidth=2, label=f'Mediana: {median_transcripts:.2f}')
plt.legend()
#plt.show()
plt.savefig('transcripts_distribution.png')
plt.clf()  # Limpar a figura atual

# Histograma de número de transcritos por função do gene
##plt.subplot(1, 2, 2)
##colors = sns.color_palette('husl', len(transcripts_per_gene_function))

#for i, (gene_function, count) in enumerate(transcripts_per_gene_function.items()):
    #sns.histplot([count], bins=range(0, transcripts_per_gene_function.max() + 10, 10), color=colors[i], label=gene_function, kde=True)
    #sns.histplot([count], bins=50, color=colors[i], label=gene_function, kde=True)  # Ajustando os bins para 50

##for gene_function, count in transcripts_per_gene_function.items():
    #sns.histplot([count], bins=50, color=colors.pop(), label=gene_function, kde=True)
    ##sns.histplot(data=data[data['Gene Function'] == gene_function], x='Gene Function', bins=range(0, transcripts_per_gene_function.max() + 10, 10), color=colors.pop(), label=gene_function, kde=True)

# OBS: esse subplot nao traz informação relevante - não plotar.
#plt.subplot(1, 2, 2)
#bx = sns.barplot(x=transcripts_per_gene_function.index, y=transcripts_per_gene_function.values, palette='husl')
#plt.title('Distribuição do número de transcritos por categoria do grupo')
#plt.xlabel('Categoria do grupo')
#plt.ylabel('Transcritos')
#plt.xticks(rotation=45)
#plt.xscale('log')
#plt.legend(title='Função do gene', bbox_to_anchor=(1.05, 1), loc='upper left')
#for index, value in enumerate(transcripts_per_gene_function.values):
#    bx.text(index, value, f'{value:,}', ha='center', va='bottom')
#plt.tight_layout()
#plt.savefig('transcripts_distribution.png')
#plt.clf()  # Limpar a figura atual

plt.figure(figsize=(10, 6))
sns.histplot(transcript_size_distribution, kde=True, bins=range(0, transcript_size_distribution.max() + 10, 10))
plt.title('Distribuição do tamanho dos transcritos')
plt.xlabel('Tamanho do transcrito')
plt.ylabel('Frequência')
plt.tight_layout()
#plt.show()
plt.savefig('transcript_size_distribution.png')
plt.clf()  # Limpar a figura atual

# 2) Quantos transcritos são lncRNAs?
lncRNA_transcripts = data[data['Transcript Function'].str.contains('lncRNA')]
num_lncRNA_transcripts = len(lncRNA_transcripts)

print(f'Número de transcritos lncRNAs: {num_lncRNA_transcripts}')

# 3) Nos genes com transcritos ncRNAs, como estão os transcritos?
#ncRNA_categories = ['ncRNA', 'lncRNA', 'protein and non-coding', 'protein and lncRNA']
ncRNA_categories = ['ncRNA', 'lncRNA']
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
#plt.show()
plt.savefig('ncRNA_transcripts_size_distribution.png')
plt.clf()  # Limpar a figura atual

# 4) Qual a distribuição dos genótipos nos genes lncRNAs? Todos os 48 genótipos estão presentes?
lncRNA_genes = lncRNA_transcripts['Gene'].unique()
genotypes_in_lncRNA_genes = data[data['Gene'].isin(lncRNA_genes)]['Genotype'].unique()

print(f'Número de genótipos nos grupos lncRNAs: {len(genotypes_in_lncRNA_genes)}')
print(f'Genótipos presentes nos grupos lncRNAs: {genotypes_in_lncRNA_genes}')

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
plt.title('Categoria do transcrito por genótipo')
plt.xlabel('Genótipo')
plt.ylabel('Frequência')
plt.tight_layout()
plt.legend(title='Categoria do transcrito')
# plt.show()
plt.savefig('transcript_function_by_genotype.png')
plt.clf()  # Limpar a figura atual

'''
# 7) Nos genes como é a distribuição dos transcritos? 
# Tem transcritos ncRNA junto com lncRNAs? São apenas lncRNAs juntos? Tem transcritos "protein" e "lncRNA" juntos?
#genes_interesse = data[data['Gene Function'].isin(ncRNA_categories)]
all_categories = ['ncRNA', 'lncRNA', 'protein and non-coding', 'protein and lncRNA', 'protein-coding']
##all_categories = ['ncRNA', 'lncRNA', 'protein and non-coding']
genes_interesse = data[data['Gene Function'].isin(all_categories)]
transcript_function_by_gene_function = genes_interesse.groupby(['Gene Function', 'Transcript Function']).size().unstack(fill_value=0)

fig, ax = plt.subplots(figsize=(10, 6))
transcript_function_by_gene_function.plot(kind='bar', stacked=True, ax=ax)
ax.set_title('Distribuição da categoria dos transcritos por categoria dos grupos')
ax.set_xlabel('Categoria do grupo')
ax.set_ylabel('Frequência dos transcritos')
plt.xticks(rotation=45)

# Adicionar anotações com os totais no topo das barras
for index, value in enumerate(transcripts_per_gene_function.values):
    ax.text(index, value, f'{value:,}', ha='center', va='bottom')

plt.tight_layout()
plt.legend(title='Categoria do transcrito')
plt.show()
#plt.savefig('transcript_function_by_gene_function.png')
plt.clf()  # Limpar a figura atual


'''

# Categorias de interesse
all_categories = ['ncRNA', 'lncRNA', 'protein and non-coding']

# Filtrando os genes de interesse
genes_interesse = data[data['Gene Function'].isin(all_categories)]

# Contando a quantidade de transcritos por Gene e Transcript Function
transcript_function_by_gene = genes_interesse.groupby(['Gene', 'Transcript Function']).size().unstack(fill_value=0)

# Agrupando por Gene Function para obter as contagens de Transcript Function dentro de cada Gene Function
gene_function_counts = genes_interesse.groupby(['Gene Function', 'Transcript Function']).size().unstack(fill_value=0)

# Plotando o gráfico de barras empilhadas
fig, ax = plt.subplots(figsize=(10, 6))
gene_function_counts.plot(kind='bar', stacked=True, ax=ax)

ax.set_title('Distribuição da categoria dos transcritos por categoria dos grupos')
ax.set_xlabel('Categoria do grupo')
ax.set_ylabel('Frequência dos transcritos')
plt.xticks(rotation=45)

# Adicionar anotações com os totais no topo das barras
for p in ax.patches:
    height = p.get_height()
    if height > 0:
        ax.annotate(f'{int(height)}', (p.get_x() + p.get_width() / 2, p.get_y() + height), ha='center', va='center', fontsize=8, color='black', xytext=(0, 5), textcoords='offset points')

plt.tight_layout()
plt.legend(title='Categoria do transcrito')
#plt.show()
plt.savefig('transcript_function_by_gene_function.png')
plt.clf()  # Limpar a figura atual
