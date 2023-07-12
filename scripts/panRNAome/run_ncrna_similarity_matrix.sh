#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./ncrna_similarity_matrix.py >> DB_clust_groups.tsv
