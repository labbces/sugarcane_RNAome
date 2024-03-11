#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./calculateClusterSimilarity.py -d ../pearson_mcl/ -t 0.5 -o Perlo
