#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 1

MATRIX=Perlo2022_counts_filters_VST_CNC_CV_above0.6.txt
MCL=Perlo2022_counts_filters_VST_CNC_CV_above0.6_mcl.txt

NODES=$(wc -l < $MATRIX)
EDGES=$(wc -l < $MCL)

echo $MATRIX has $NODES nodes and $EDGES edges.
