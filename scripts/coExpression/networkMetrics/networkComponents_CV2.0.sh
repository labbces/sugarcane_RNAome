#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 1

MATRIX=Correr2020_counts_filters_VST_CNC_CV_above_2.txt
MCL=Correr2020_counts_filters_VST_CNC_CV_above2_mcl.txt

NODES=$(wc -l < $MATRIX)
EDGES=$(wc -l < $MCL)

echo $MATRIX has $NODES nodes and $EDGES edges.
