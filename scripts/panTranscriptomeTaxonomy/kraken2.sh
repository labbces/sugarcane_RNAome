#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 20
#$ -l h=figsrv

INFILE=`ls -1 48_transcriptomes.fasta`
OUTBASE=`basename $INFILE`
OUTPUT=${OUTBASE/.fasta}.out.kraken
REPORT=${OUTBASE/.fasta}.report.kraken

module load Kraken2/2.1.2

kraken2 --threads $NSLOTS --confidence 0.05 --report-zero-counts --report $REPORT --output $OUTPUT --db /Storage/data1/felipe.peres/Kraken2/krakendb $INFILE
