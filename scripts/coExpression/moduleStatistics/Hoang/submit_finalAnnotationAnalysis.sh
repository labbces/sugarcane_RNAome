#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./finalAnnotationAnalysis.py -annotation guiltByAssociationFunctions.tsv -obo ../enrichedGOFrequency/go-basic.obo -out_ncRNA ncRNAs_GOFrequencies.tsv -out_lncRNA lncRNAs_GOFrequencies.tsv -plot lncRNAs_sizeDistribution.png
