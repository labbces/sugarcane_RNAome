#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 10

module load miniconda3
conda activate networkx

network=`ls -1 *_mcl.txt`
out=${network/.txt/_degree.tsv}

# columns = source node, target node, weight
/usr/bin/time -v ./parallelBetweennessCentrality.py --edge_list $network --output $out --threads $NSLOTS --columns 0 1 2
