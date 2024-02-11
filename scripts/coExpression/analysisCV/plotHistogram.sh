#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -t 1
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./HistogramByClassificationAndFunction.py -cv Hoang2017_counts_filters_VST_CV_classified.txt -oc Hoang2017_filtered_CV_Classification.png -of Hoang2017_filtered_CV_Function.png -b 50 -t 1 1.5 1.8 2.0
