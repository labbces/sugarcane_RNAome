#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -t 1
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./classifyGenesfromCVfile.py -c panTranscriptome_panRNAomeClassificationTable_hyphen_Class.tsv -cv Hoang2017_counts_filters_VST_CV.txt -o Hoang2017_counts_filters_VST_CV_classified.txt
