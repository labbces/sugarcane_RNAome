#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -t 1
#$ -pe smp 1

module load miniconda3
conda activate plotly

module load R/4.0.0

/usr/bin/time -v Rscript speciesOfOrigin_filter1.8CPM.R
