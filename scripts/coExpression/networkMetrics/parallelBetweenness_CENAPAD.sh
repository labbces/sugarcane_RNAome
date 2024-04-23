#!/bin/bash

#PBS -N parallelBetweenness
#PBS -q memlong
#PBS -l nodes=1:ppn=128
#PBS -e parallelBetweenness.sh.err
#PBS -o parallelBetweenness.sh.out

cd $PBS_O_WORKDIR

source /home/lovelace/proj/proj832/fvperes/miniconda3/etc/profile.d/conda.sh
conda activate networkx

network=`ls -1 *_mcl.txt`
out=${network/.txt/_betweenness.tsv}

# columns = source node, target node, weight
/usr/bin/time -v ./parallelBetweennessCentrality.py --edge_list $network --output $out --threads 128 --columns 0 1 2
