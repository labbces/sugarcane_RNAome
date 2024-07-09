#!/usr/bin/env python

import pandas as pd

annotation_file = 'updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv'
annotation_df = pd.read_csv(annotation_file, sep='\t', header=None, low_memory=False)
annotation_df.columns = ['panRNAome_category', 'gene_name', 'transcript_name', 'gene_function', 'transcript_length', 'transcript_function', 'Transcript Rfam family', 'Transcript GO', 'Gene GO']

origin_file = 'speciesOfOriginPanTranscriptome_FPKM.tsv'
origin_df = pd.read_csv(origin_file, sep='\t', header=0)

# mapear cada transcrito ao genotipo de origem
transcript_to_origin = dict(zip(origin_df['Transcript'], origin_df['Origin']))

# nova coluna 'Transcript Origin' ao arquivo de anotacao
annotation_df['Transcript Origin'] = annotation_df['transcript_name'].map(transcript_to_origin)

# agrupar a origem todos os transcritos pelo grupo - separados por virgula
group_to_origin = annotation_df.groupby('gene_name')['Transcript Origin'].apply(lambda x: ','.join(x.dropna().unique())).to_dict()

# nova coluna 'Gene Origin' ao arquivo de anotacao
annotation_df['Gene Origin'] = annotation_df['gene_name'].map(group_to_origin)

# Convertendo transcript_length para inteiro (removendo ".0" se for inteiro)
annotation_df['transcript_length'] = annotation_df['transcript_length'].astype('Int64')  # Utiliza tipo Int64 para manter NA e inteiros

output_file = 'updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO_with_Origin.tsv'
annotation_df.to_csv(output_file, sep='\t', index=False, header=None)
