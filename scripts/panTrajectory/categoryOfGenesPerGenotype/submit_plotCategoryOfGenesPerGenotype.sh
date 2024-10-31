#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./plotCategoryOfGenesPerGenotype.py -a updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv -c Exclusive
/usr/bin/time -v ./plotCategoryOfGenesPerGenotype.py -a updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv -c Soft-core
/usr/bin/time -v ./plotCategoryOfGenesPerGenotype.py -a updated_panTranscriptome_panRNAome_GeneFunction_Length_with_Rfam_with_GO.tsv -c Accessory
