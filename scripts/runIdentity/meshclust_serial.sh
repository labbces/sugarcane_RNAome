#!/bin/bash

#PBS -N MESHCLUST
#PBS -q serial
#PBS -l mem=10G
#PBS -e meshclust.sh.err
#PBS -o meshclust.sh.out

cd $PBS_O_WORKDIR

source /home/lovelace/proj/proj832/fvperes/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake -p -s Snakefile --cores 16
