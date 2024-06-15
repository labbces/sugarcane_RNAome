#!/usr/bin/env python

import argparse
import csv

parser = argparse.ArgumentParser(prog='mergeLengthWithClassificationForAllTranscripts.py', description='transfer length for all transcripts', add_help=True)
parser.add_argument('-length', dest='length_file', metavar='<48_transcriptomes_length.tsv>',help="transcripts length file", required=True)
parser.add_argument('-classification', dest='classified_file', metavar='<panTranscriptome_panRNAomeClassificationTable_hyphen_Class_length.tsv>',help="panRNAome classification with length", required=True)

args = parser.parse_args()
classified_file = args.classified_file
length_file = args.length_file
output_file = "updated_panTranscriptome_panRNAome_Length.tsv"

length_data = {}

with open(length_file, 'r') as infile:
    reader = csv.reader(infile, delimiter='\t')
    next(reader)
    for row in reader:
        sequence_id = row[0]
        length = row[1]
        length_data[sequence_id] = length

with open(classified_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.reader(infile, delimiter='\t')
    writer = csv.writer(outfile, delimiter='\t')
    
    header = next(reader)
    writer.writerow(header)
    
    for row in reader:
        sequence_id = row[2]
        length = row[4]
        type_column = row[3]
        
        if length == 'NA':
            if sequence_id in length_data:
                row[4] = length_data[sequence_id]
            if type_column == 'protein-coding':
                row[5] = 'protein-coding'
        
        writer.writerow(row)
