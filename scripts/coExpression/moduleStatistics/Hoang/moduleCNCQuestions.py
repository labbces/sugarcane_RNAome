#!/usr/bin/env python

import pandas as pd

df = pd.read_csv('moduleSummary.tsv', sep='\t')  

# agrupar as colunas 'Genes_coding' e 'Genes_protein_and_non_coding'
df['Genes_coding'] = df['Genes_coding'] + df['Genes_protein_and_non_coding']

# criar um novo df com as novas colunas
df_novo = df[['Module', 'Genes', 'Genes_coding', 'Genes_ncRNA', 'Genes_lncRNA']]

# Quantos modulos tem CNC?
# (modulos com pelo menos 1 gene nas colunas 'Genes_coding', 'Genes_ncRNA', e 'Genes_lncRNA')
modulos_cnc = df_novo[
    (df_novo['Genes_coding'] > 0) & 
    ((df_novo['Genes_ncRNA'] > 0) | (df_novo['Genes_lncRNA'] > 0))
]

qtd_modulos_cnc = modulos_cnc.shape[0]

# Quantos módulos CNC possuem Genes_lncRNA?
modulos_cnc_com_lncRNA = modulos_cnc[modulos_cnc['Genes_lncRNA'] > 0]
qtd_modulos_cnc_com_lncRNA = modulos_cnc_com_lncRNA.shape[0]

modulos_cnc_com_lncRNA.to_csv('moduleCNCwithlncRNAs.tsv', index=False, sep='\t')

with open('moduleCNCAnswer.txt', 'w') as f:
	f.write(f"1) Quantos módulos totais? {df_novo.shape[0]}\n")
	f.write(f"2) Quantos módulos são CNC? {qtd_modulos_cnc}\n")
	f.write(f"3) Quantos módulos CNC possuem Genes_lncRNA? {qtd_modulos_cnc_com_lncRNA}\n")
