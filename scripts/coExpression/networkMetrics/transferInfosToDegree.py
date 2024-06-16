#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(prog='transferInfosToDegree.py', description='transfer panRNAome classification to degree file', add_help=True)
parser.add_argument('-classification', dest='input_classification_table', metavar='<updated_panTranscriptome_panRNAome_Length.tsv>',help="updated panRNAome classification", required=True)
parser.add_argument('-degree', dest='input_degree_file', metavar='<degree.tsv>',help="file with genes degree", required=True)

args = parser.parse_args()
input_classification_table = args.input_classification_table
input_degree_file = args.input_degree_file
output_file = input_degree_file[:-4] + "_updated.tsv"

info_data = {}

#Soft-core       OG0000000       B1_k25_TRINITY_DN12555_c1_g1_i10        protein-coding  309     protein-coding
with open(input_classification_table, 'r') as infile:
    for line in infile:
        parts = line.strip().split('\t')
        panRNAome = parts[0]
        gene = parts[1]
        transcript = parts[2]
        gene_function = parts[3]
        transcript_length = parts[4]
        transcript_function = parts[5]
        info_data[gene] = (panRNAome, transcript, gene_function, transcript_length, transcript_function)

with open(input_degree_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        row = line.strip().split('\t')
        gene = row[0]
        
        if gene in info_data:
            combined_row = line.strip() + '\t' + '\t'.join(info_data[gene])
            outfile.write(combined_row + '\n')
