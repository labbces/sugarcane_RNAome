#!/usr/bin/env python

import argparse
import csv

parser = argparse.ArgumentParser(prog='mergeLengthWithClassificationForAllTranscripts.py', description='transfer length for all transcripts', add_help=True)
parser.add_argument('-length', dest='length_file', metavar='<48_transcriptomes_length.tsv>',help="transcripts length file", required=True)
parser.add_argument('-classification', dest='classified_file', metavar='<panTranscriptome_panRNAomeClassificationTable_hyphen_Class_length.tsv>',help="panRNAome classification with length", required=True)
parser.add_argument('-protein_coding', dest='protein_coding_file', metavar='<proteinCoding_RNA.txt>', help="file with protein-coding transcripts", required=True)

args = parser.parse_args()
classified_file = args.classified_file
length_file = args.length_file
protein_coding_file = args.protein_coding_file
output_file = "updated_panTranscriptome_panRNAome_Length.tsv"

length_data = {}
protein_coding_transcripts = set()

# Load the length data
with open(length_file, 'r') as infile:
    reader = csv.reader(infile, delimiter='\t')
    next(reader)  # Skip header
    for row in reader:
        sequence_id = row[0]
        length = row[1]
        length_data[sequence_id] = length

# Load the protein-coding transcripts
with open(protein_coding_file, 'r') as infile:
    for line in infile:
        transcript_id = line.strip()
        protein_coding_transcripts.add(transcript_id)

# Process the classified file and write to the output file
with open(classified_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.reader(infile, delimiter='\t')
    writer = csv.writer(outfile, delimiter='\t')
    
    header = next(reader)
    writer.writerow(header)
    
    for row in reader:
        sequence_id = row[2]
        current_length = row[4]
        current_function = row[5]
        
        # Update length if it is 'NA'
        if current_length == 'NA' and sequence_id in length_data:
            row[4] = length_data[sequence_id]
        
        # Update transcript function
        if sequence_id in protein_coding_transcripts:
            if current_function == 'ncRNA':
                row[5] = 'protein and non-coding'
            elif current_function == 'lncRNA':
                row[5] = 'protein and lncRNA'
            else:
                row[5] = 'protein-coding'
        
        writer.writerow(row)

