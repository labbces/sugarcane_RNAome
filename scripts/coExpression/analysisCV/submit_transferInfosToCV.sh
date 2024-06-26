#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 1

module load miniconda3
/usr/bin/time -v ./transferInfosToCV.py -classification updated_panTranscriptome_panRNAome_GeneFunction_Length.tsv -cv Correr2020_counts_filters_VST_CV.txt
