#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(prog='finalAnnotationGuiltByAssociation.py', description='generate lnc/ncRNAs annotation using guilt by association', add_help=True)
parser.add_argument('-genes', dest='genes_go_file', metavar='genes_ncRNA_com_enrichedGO_final.tsv', help='lnc/ncRNAs and associated GOs', type=str, required=True)
parser.add_argument('-annotation', dest='annotation_file', metavar='updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv', help='Annotation table', type=str, required=True)
parser.add_argument('-out', dest='output_file', metavar='guiltByAssociationFunctions.tsv', help='Output table with GO associated to lnc/ncRNAs', type=str, required=True)
args = parser.parse_args()

#genes_go_file = 'genes_ncRNA_com_enrichedGO_final.tsv'
#annotation_file = '../updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv'
#output_file = 'guiltByAssociationFunctions.tsv'

genes_go_file = args.genes_go_file
annotation_file = args.annotation_file
output_file = args.output_file

def combine_go_annotations(genes_go_file, annotation_file, output_file):
    genes_go_df = pd.read_csv(genes_go_file, sep='\t')
    
    annotation_df = pd.read_csv(annotation_file, sep='\t', index_col=False)

    filtered_annotation_df = annotation_df[annotation_df['Gene'].isin(genes_go_df['Gene'])]

    combined_df = pd.merge(filtered_annotation_df, genes_go_df[['Gene', 'GO (Guilt by association)']],
                           on='Gene', how='inner')

    combined_df.to_csv(output_file, sep='\t', index=False)
    print(f"Guilt by association table saved: {output_file}")

combine_go_annotations(genes_go_file, annotation_file, output_file)
