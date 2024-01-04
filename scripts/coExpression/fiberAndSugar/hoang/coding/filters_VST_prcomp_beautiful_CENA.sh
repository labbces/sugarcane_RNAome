#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -t 1
#$ -pe smp 1

module load R/4.0.0

/usr/bin/time -v Rscript redes_hoang_filters_VST_prcomp_coding_CENA.R
