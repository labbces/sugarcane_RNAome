#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./getGOTermsUsedToAnnotateEnrichedModules.py -modules ../enrichedGOFrequency/moduleCNCwithlncRNAswithEnrichedGO.tsv -genes ../Correr2020_counts_filters_VST_topCV_mcl_formated_cliques.csv -annotation ../updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv -out GOTermsUsedInModuleCNCwithlncRNAswithEnrichedGO.tsv
