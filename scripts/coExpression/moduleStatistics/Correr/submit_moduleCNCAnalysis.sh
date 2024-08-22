#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./moduleCNCAnalysis.py -m Correr2020_counts_filters_VST_topCV_mcl_formated_cliques.csv -f Correr2020_counts_filters_VST_topCV_mcl_formated_cliques_filtered.csv -a updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv -o moduleSummary.tsv
