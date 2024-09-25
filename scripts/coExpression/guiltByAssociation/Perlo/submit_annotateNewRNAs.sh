#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./annotateNewRNAs.py -m ../Perlo2022_counts_filters_VST_topCV_mcl_formated_cliques.csv -go ../enrichedGOFrequency -a ../updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv -o genes_ncRNA_com_GO_final.tsv
