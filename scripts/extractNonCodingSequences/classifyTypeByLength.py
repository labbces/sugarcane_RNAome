#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(prog='classifyTypeByLength.py', description='classify ncRNA by length', add_help=True)
parser.add_argument('-i', dest='input_file', metavar='<putative_ncRNAs_length.tsv>',help="ncRNAs length file", required=True)

args = parser.parse_args()
input_file = args.input_file
output_file = input_file[:-4] + "_classified.tsv"

def classify_sequence(length):
    return 'lncRNA' if length >= 500 else 'ncRNA'

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    
    header = infile.readline().strip()
    outfile.write(header + '\tType\n')
    
    for line in infile:
        parts = line.strip().split('\t')
        sequence_id = parts[0]
        length = int(parts[1])
        classification = classify_sequence(length)
        outfile.write(f"{sequence_id}\t{length}\t{classification}\n")
