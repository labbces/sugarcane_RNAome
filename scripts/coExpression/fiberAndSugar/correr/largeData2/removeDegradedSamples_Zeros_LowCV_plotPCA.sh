#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -t 1
#$ -pe smp 1

module load R/4.0.0

/usr/bin/time -v Rscript redes_correr_tximport_calcZeros_VST_plotPCA_CENA.R
