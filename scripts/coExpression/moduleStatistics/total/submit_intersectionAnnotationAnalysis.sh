#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./finalAnnotationAnalysis.py -annotation intersection_guiltByAssociationFunctions.tsv -obo go-basic.obo -out_ncRNA intersection_ncRNAs_GOFrequencies.tsv -out_lncRNA intersection_lncRNAs_GOFrequencies.tsv -plot intersection_lncRNAs_sizeDistribution.png
