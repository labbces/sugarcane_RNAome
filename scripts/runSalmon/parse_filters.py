#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(prog='parse_filters.py', description='Parse CoNekT Bioenergy snakemake pipeline filters', add_help=True)
parser.add_argument('--accessions', dest='accessions', metavar='<accesions csv file>', help='file containing SRA IDs', required=True)
parser.add_argument('--paired_srrlist', dest='paired_srrlist', metavar='<paired srrlist csv file>', help='file containing only paired accessions', required=True)
parser.add_argument('--filter_stats', dest='filter_stats', metavar='<filter stats file>', help='file containing only filtered accessions', required=True)

args = parser.parse_args()
accessions = args.accessions
paired_srrlist = args.paired_srrlist
filter_stats = args.filter_stats

# trocando virgula por \n no acessions (dentro do shell snakemake)
sample_column = ['Sample']
download = pd.read_csv(accessions, header=None, names=sample_column)

# trocando virgula por \n no paired_srrlist (dentro do shell snakemake)
filter1 = pd.read_csv(paired_srrlist, header=None, names=sample_column)

# trocando tab por virgula no filter_stats
filter_columns = ['Sample', 'Mapping Rate', 'Perc Low TPM']
filter2 = pd.read_csv(filter_stats, header=0, names=filter_columns)

# concatenando samples
merged = pd.merge(download, filter2, on='Sample', how='left')
merged['Downloaded'] = merged['Sample'].isin(filter2['Sample'])
merged['Pass Filter'] = merged['Sample'].isin(filter1['Sample'])

merged.to_csv('preliminar_report.tsv', index=False)
print(merged)
