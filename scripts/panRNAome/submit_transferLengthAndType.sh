#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./transferLengthAndType.py -length putative_ncRNAs_length_classified.tsv -classification panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv
