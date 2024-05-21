#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 1

MATRIX=Correr2020_counts_filters_VST_CNC_CV_above_1.5.txt
MCL=Correr2020_counts_filters_VST_CNC_CV_above_1.5_mcl_combined.txt

NODES=$(wc -l < $MATRIX)
EDGES=$(wc -l < $MCL)

echo $MATRIX has $NODES nodes and $EDGES edges.
