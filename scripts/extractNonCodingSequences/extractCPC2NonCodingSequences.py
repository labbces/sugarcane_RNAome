#!/usr/bin/env python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(prog='extractCPC2NonCodingSequences.py', description='extract sequences labelled as noncoding by CPC2', add_help=True)
parser.add_argument('-cpcout', dest='cpcout', metavar='<CPC2 output (delimited by tab)>',help="CPC2 result in table format", required=True)
parser.add_argument('-transcriptome', dest='transcriptome_file', metavar='<transcriptome.fasta>',help="transcriptome file", required=True)

args = parser.parse_args()
cpcout = args.cpcout
transcriptome_file = args.transcriptome_file

output_file = transcriptome_file[:-5] + "cpc_ncrnas.fa"

#CPCout = CPC2 result in table format (delimited by tab):
#ID peptide_length Fickett_score isoelectric_point ORF_integrity coding_probability coding_label

noncoding_list = []
with open(cpcout, "r") as f:
	for line in f:
		ID = line.split()[0]
		coding_label = line.split()[7]
		if coding_label == "noncoding":
			noncoding_list.append(ID)
#print(noncoding_list)

fin = open(transcriptome_file, "r")
fout = open(output_file, "w")

for record in SeqIO.parse(fin, "fasta"):
	for item in noncoding_list:
		if item.strip() == record.id:
			SeqIO.write(record, fout, "fasta")
