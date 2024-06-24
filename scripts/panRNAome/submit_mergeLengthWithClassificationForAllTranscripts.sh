#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./mergeLengthWithClassificationForAllTranscripts.py -length 48_transcriptomes_length.tsv -classification panTranscriptome_panRNAomeClassificationTable_hyphen_Class_length.tsv -protein_coding proteinCoding_RNA.txt
