#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 1

OrthoFile=`ls DB_clust_groups_withoutLastColumn.tsv`

module load miniconda3

/usr/bin/time -v ./extract_pan_transcriptome_data.py -i $OrthoFile
