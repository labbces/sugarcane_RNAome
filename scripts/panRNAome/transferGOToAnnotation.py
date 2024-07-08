#!/usr/bin/env python

import pandas as pd

annotation_file = 'updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam.tsv'
annotation_df = pd.read_csv(annotation_file, sep='\t', header=None, low_memory=False)
annotation_df.columns = ['panRNAome_category', 'gene_name', 'transcript_name', 'gene_function', 'transcript_length', 'transcript_function', 'Transcript Rfam family']

GO_file = 'GO_universe_annotation_list'

# alguns headers foram duplicados pois o arquivo com GO_universe é resultado do GO_universe de cada um dos 48 genotipos
# a linha dos seguintes headers estao duplicadas (remover -1 de cada linha, pois panda comeca contar do 0)
header_lines_to_skip = [3161324,6305241,9743264,13115429,16753257,20697789,24112970,27329728,30968518,34592847,37780558,39570638,41175849,43956902,46868010,49557974,51191518,53320534,55532197,57165272,59157312,60599337,62828437,64828109,69997279,73066428,78120993,80627945,82588607,85252039,87930459,89935273,92151508,94495039,96926290,99252945,101218452,103594818,106482027,108966306,113597443,118480642,121307191,123314414,126310238,128484350,130971502]

GO_df = pd.read_csv(GO_file, sep='\t', header=0, skiprows=header_lines_to_skip, dtype={
	'qpid': str, 
	'ontology': str, 
	'goid': str, 
	'desc': str, 
	'ARGOT_score': float, 
	'ARGOT_PPV': float, 
	'ARGOT_rank': int, 
	'goclasscount': int
}) 

# filtrar ARGOT_PPV > 0.6 e ontology == 'BP'
filtered_GO_df = GO_df[(GO_df['ARGOT_PPV'] > 0.6) & (GO_df['ontology'] == 'BP')]

# adicionar "GO:" antes dos identificadores goid
filtered_GO_df['goid'] = 'GO:' + filtered_GO_df['goid']

# mapear cada transcrito ao Gene Ontology
transcript_to_go = dict(zip(filtered_GO_df['qpid'], filtered_GO_df['goid']))

# nova coluna 'Transcript GO (BP)' ao arquivo de anotacao
annotation_df['Transcript GO (BP)'] = annotation_df['transcript_name'].map(transcript_to_go)

# filtrar arquivo de anotacao para adicionar a coluna Gene GO (BP) apenas para genes protein-coding
protein_df = annotation_df[annotation_df['gene_function'] == 'protein-coding']

# agrupar os GO de todos os transcritos por gene - separados por virgula
#gene_to_go = annotation_df.groupby('gene_name')['Transcript GO (BP)'].apply(lambda x: ','.join(x.dropna().unique())).to_dict()
gene_to_go = protein_df.groupby('gene_name')['Transcript GO (BP)'].apply(lambda x: ','.join(x.dropna().unique())).to_dict()

# nova coluna 'Gene GO (BP)' ao arquivo de anotacao
annotation_df['Gene GO (BP)'] = annotation_df['gene_name'].map(gene_to_go)

# Convertendo transcript_length para inteiro (removendo ".0" se for inteiro)
annotation_df['transcript_length'] = annotation_df['transcript_length'].astype('Int64')  # Utiliza tipo Int64 para manter NA e inteiros

output_file = 'updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv'
annotation_df.to_csv(output_file, sep='\t', index=False, header=None)
