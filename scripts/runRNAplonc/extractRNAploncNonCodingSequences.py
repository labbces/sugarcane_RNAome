#!/usr/bin/env python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(prog='extractRNAploncNonCodingSequences.py', description='extract sequences labelled as noncoding by RNAplonc', add_help=True)
parser.add_argument('-rnaplonc', dest='rnaploncout', metavar='<RNAplonc output (delimited by tab)>',help="RNAplonc result in table format", required=True)
parser.add_argument('-transcriptome', dest='transcriptome_file', metavar='<transcriptome.fasta>',help="transcriptome file", required=True)

args = parser.parse_args()
rnaploncout = args.rnaploncout
transcriptome_file = args.transcriptome_file

output_file = transcriptome_file[:-3] + "_rnaplonc_ncrnas.fa"

#RNAploncout = RNAplonc result in table format (delimited by tab):
#example file

noncoding_list = []
with open(rnaploncout, "r") as f:
	for line in f:
		ID = line.split()[0]
		noncoding_list.append(ID)
#print(noncoding_list)

fin = open(transcriptome_file, "r")
fout = open(output_file, "w")

for record in SeqIO.parse(fin, "fasta"):
	for item in noncoding_list:
		if item.strip() == record.id:
			SeqIO.write(record, fout, "fasta")
