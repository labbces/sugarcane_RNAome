#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 1

module load miniconda3
/usr/bin/time -v ./calculateDegree.py -i Correr2020_counts_filters_VST_CNC_CV_above2_mcl.txt
