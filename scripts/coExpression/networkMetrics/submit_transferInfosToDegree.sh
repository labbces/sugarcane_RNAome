#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 1

module load miniconda3
/usr/bin/time -v ./transferInfosToDegree.py -classification updated_panTranscriptome_panRNAome_Length.tsv -degree Correr2020_counts_filters_VST_CNC_CV_above_1.5_mcl_combined_degree_sorted.tsv

/usr/bin/time -v ./transferInfosToDegree.py -classification updated_panTranscriptome_panRNAome_Length.tsv -degree Correr2020_counts_filters_VST_CNC_CV_above2_mcl_degree_sorted.tsv
