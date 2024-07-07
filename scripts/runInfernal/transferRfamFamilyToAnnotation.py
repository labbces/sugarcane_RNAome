#!/usr/bin/env python

import pandas as pd

annotation_file = 'updated_panTranscriptome_panRNAome_GeneFunction_Length.tsv'
annotation_df = pd.read_csv(annotation_file, sep='\t', header=None)
annotation_df.columns = ['panRNAome_category', 'gene_name', 'transcript_name', 'gene_function', 'transcript_length', 'transcript_function']

rfam_file = 'putative_ncRNA_consensus.complete.deoverlapped.cmscan.tab.tblout'
#cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26, 27]
rfam_columns = ['idx', 'target_name', 'target_accession', 'query_name', 'query_accession', 'clan', 'name', 'mdl', 'mdl_from', 'mdl_to', 'seq_from', 'seq_to', 'strand', 'trunc', 'pass', 'gc', 'bias', 'score', 'E-value', 'inc', 'olp', 'anyidx', 'afrct1', 'afrct2', 'winidx', 'wfrct1', 'wfrct2']

# obs: ignorar coluna com descrição da Rfam family - coluna 28
rfam_df = pd.read_csv(rfam_file, sep='\t', header=None, usecols=range(27), names=rfam_columns)

# mapear cada transcrito à sua família Rfam
transcript_to_rfam = dict(zip(rfam_df['query_name'], rfam_df['target_accession']))

# nova coluna 'Transcript Rfam family' ao arquivo de anotação
annotation_df['Transcript Rfam family'] = annotation_df['transcript_name'].map(transcript_to_rfam)

# agrupar as Rfam family de todos os transcritos por gene - separados por virgula e.g: RF04085,RF04087,RF04089
#gene_to_rfam = annotation_df.groupby('gene_name')['Transcript Rfam family'].apply(lambda x: ','.join(x.dropna().unique())).to_dict()

# nova coluna 'Gene Rfam family' ao arquivo de anotação
#annotation_df['Gene Rfam family'] = annotation_df['gene_name'].map(gene_to_rfam)

# Convertendo transcript_length para inteiro (removendo ".0" se for inteiro)
annotation_df['transcript_length'] = annotation_df['transcript_length'].astype('Int64')  # Utiliza tipo Int64 para manter NA e inteiros

output_file = 'updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam.tsv'
annotation_df.to_csv(output_file, sep='\t', index=False, header=None)
