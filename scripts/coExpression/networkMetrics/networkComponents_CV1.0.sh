#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 1

MATRIX=Hoang2017_counts_filters_VST_CNC_CV_above_1.txt
MCL=Hoang2017_counts_filters_VST_CNC_CV_above_1_mcl.txt

NODES=$(wc -l < $MATRIX)
EDGES=$(wc -l < $MCL)

echo $MATRIX has $NODES nodes and $EDGES edges.
