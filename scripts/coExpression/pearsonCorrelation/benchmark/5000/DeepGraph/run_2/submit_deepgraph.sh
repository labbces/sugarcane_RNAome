#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 6

source /Storage/progs/miniconda3/etc/profile.d/conda.sh
conda activate deepgraph

INFILE=`ls -1 *0.txt` 

# Calculate pearson correlation 
# n_process = 6
/usr/bin/time -v ./deepgraph_pearsonCorrelation.py -i $INFILE

# Parse output and generate table of correlations
/usr/bin/time -v ./deepgraph_parseOutput.py -i $INFILE

