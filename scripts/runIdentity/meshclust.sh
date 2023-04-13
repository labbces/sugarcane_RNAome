#!/bin/bash

#PBS -N MESHCLUST
#PBS -q par16
#PBS -l nodes=1:ppn=16
#PBS -e meshclust.sh.err
#PBS -o meshclust.sh.out

cd $PBS_O_WORKDIR

source /home/lovelace/proj/proj832/fvperes/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake -p -s Snakefile --cores 16
