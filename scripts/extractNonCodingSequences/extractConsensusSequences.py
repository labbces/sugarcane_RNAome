#!/usr/bin/env python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(prog='extractConsensusSequences.py', description='extract consensus sequences by list from fasta', add_help=True)
parser.add_argument('-transcriptome', dest='transcriptome_file', metavar='<transcriptome.fasta>',help="transcriptome file", required=True)
parser.add_argument('-list', dest='list_file', metavar='<sequences_to_extract.txt>', help="file containing sequences to extract", required=True)

args = parser.parse_args()
transcriptome_file = args.transcriptome_file
list_file = args.list_file

output_file = "putative_ncRNA_consensus.fa"

# Read the file containing sequences to extract and create a set to speed up membership verification
sequences_to_extract_set = set()
with open(list_file, "r") as list_f:
    for line in list_f:
        sequence = line.strip()
        sequences_to_extract_set.add(sequence)

with open(transcriptome_file, "r") as fin, open(output_file, "w") as fout:
    for record in SeqIO.parse(fin, "fasta"):
        if record.id in sequences_to_extract_set:
            SeqIO.write(record, fout, "fasta")
