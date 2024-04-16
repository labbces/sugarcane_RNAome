#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 1

network=`ls -1 *_mcl.txt`

# columns = source node, target node, weight
/usr/bin/time -v ./parallelBetweennessCentrality.py --edge_list $network --output test_betweenness.csv --threads 2 --columns 0 1 2
