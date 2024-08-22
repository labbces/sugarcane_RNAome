#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./moduleCNCAnalysis.py -m Hoang2017_counts_filters_VST_topCV_mcl_formated_cliques.csv -f Hoang2017_counts_filters_VST_topCV_mcl_formated_cliques_filtered.csv -a updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv -o moduleSummary.tsv
