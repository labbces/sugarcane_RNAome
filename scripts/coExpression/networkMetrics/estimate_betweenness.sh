#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 1

network=`ls -1 *_mcl.txt`

module load R/4.0.0

Rscript network_analysis_estimate_betweenness.r $network
