#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 1

INFILE=DB_clust_groups_withoutLastColumn.tsv

module load Python/3.7.2
/usr/bin/time -v python3.7 extractPanTranscriptomeGroups.py -i ${INFILE} -p 0.8
