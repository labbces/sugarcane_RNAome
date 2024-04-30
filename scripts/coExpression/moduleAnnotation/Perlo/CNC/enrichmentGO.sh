#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -t 1
#$ -pe smp 1

module load R/4.0.0

/usr/bin/time -v Rscript enrichmentGO.R clusterSize.txt GO_annotations_BP_PPV0.6.tsv Perlo2022_counts_filters_VST_topCV_mcl_formated_cliques.csv 
