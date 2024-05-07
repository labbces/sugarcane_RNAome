#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 1

INFILE=DB_clust_Orthogroups.tsv

module load Python/3.7.2
/usr/bin/time -v python3.7 computeOrthogroupStats_0.8.py --orthogroupsFile ${INFILE} --numberSpecies 48 --suffixOut 0.8
