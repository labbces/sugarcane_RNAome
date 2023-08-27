#!/usr/bin/env python

import gzip
from Bio import SeqIO
import argparse
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
                    help="provide the name of the sequence input file")
parser.add_argument("format", type=str,
                    help="provide the format of the sequences in your input file")
parser.add_argument("prefixOut", type=str,
                    help="provide a prefix to create output files")

args=parser.parse_args()

pgcList=[]
lenList=[]
totalReads=0
totalLength=0
print(args.infile)
with open(args.infile, "rt") as handle:
    for record in SeqIO.parse(handle, args.format):
#        print(record.id)
        pgc=(record.seq.upper().count('G') + record.seq.upper().count('C')) / len(record.seq)
        pgcList.append(pgc)
        lenList.append(len(record.seq))
        totalLength=totalLength+len(record.seq)
        totalReads=totalReads+1

fig, (ax1, ax2,ax3) = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))
ax1.hist(lenList,bins=1000,histtype='bar')
ax1.set_title('Read lengths')
ax1.set_xlabel('Length (bp)')
ax1.tick_params(axis='x', labelsize=8)
ax1.tick_params(axis='y', labelsize=8)

ax2.hist(lenList,bins=1000,histtype='bar')
ax2.set_title('Read lengths (log)')
ax2.set_xlabel('Length (log bp)')
ax2.set_xscale('log')
ax2.yaxis.set_ticklabels([])
ax2.tick_params(axis='x', labelsize=8)

ax3.hist(pgcList,bins=100,histtype='bar')
ax3.set_title('GC content')
ax3.set_xlabel('GC content')
ax3.yaxis.tick_right()
ax3.tick_params(axis='x', labelsize=8)
ax3.tick_params(axis='y', labelsize=8)

histogramsOutSVG=args.prefixOut+'_histograms.svg'
histogramsOutPNG=args.prefixOut+'_histograms.png'

plt.savefig(histogramsOutSVG, format="svg")
plt.savefig(histogramsOutPNG, format="png")

print(f'Total Length: {totalLength}; TotalReads={totalReads}')
