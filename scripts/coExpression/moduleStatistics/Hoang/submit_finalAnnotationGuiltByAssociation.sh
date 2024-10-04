#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./finalAnnotationGuiltByAssociation.py -genes genes_ncRNA_com_enrichedGO_final.tsv -annotation ../updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv -out guiltByAssociationFunctions.tsv 
