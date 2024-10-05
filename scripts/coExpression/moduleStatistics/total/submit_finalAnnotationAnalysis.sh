#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./finalAnnotationAnalysis.py -annotation total_guiltByAssociationFunctions.tsv -obo go-basic.obo -out_ncRNA total_ncRNAs_GOFrequencies.tsv -out_lncRNA total_lncRNAs_GOFrequencies.tsv -plot total_lncRNAs_sizeDistribution.png
