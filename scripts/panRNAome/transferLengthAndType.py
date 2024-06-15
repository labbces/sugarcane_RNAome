#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(prog='transferLengthAndType.py', description='transfer length and type to panRNAome classification', add_help=True)
parser.add_argument('-length', dest='input_classified_file', metavar='<putative_ncRNAs_length_classified.tsv>',help="ncRNAs length file", required=True)
parser.add_argument('-classification', dest='input_classification_table', metavar='<panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv>',help="panRNAome classification", required=True)

args = parser.parse_args()
input_classified_file = args.input_classified_file
input_classification_table = args.input_classification_table
output_file = input_classification_table[:-4] + "_length.tsv"

classified_data = {}

with open(input_classified_file, 'r') as infile:
    infile.readline()
    for line in infile:
        parts = line.strip().split('\t')
        sequence_id = parts[0]
        length = parts[1]
        seq_type = parts[2]
        classified_data[sequence_id] = (length, seq_type)

with open(input_classification_table, 'r') as infile, open(output_file, 'w') as outfile:
    #header = infile.readline().strip() + '\tLength\tType\n'
    #outfile.write(header)
    
    for line in infile:
        parts = line.strip().split('\t')
        sequence_id = parts[2]
        
        if sequence_id in classified_data:
            length, seq_type = classified_data[sequence_id]
        else:
            length, seq_type = 'NA', 'NA'
        
        new_line = line.strip() + f'\t{length}\t{seq_type}\n'
        outfile.write(new_line)
