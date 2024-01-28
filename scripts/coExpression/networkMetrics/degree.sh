#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 1

module load R/4.0.0
Rscript network_analysis_degree.r
