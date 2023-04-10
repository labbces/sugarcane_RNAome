#!/bin/bash

#PBS -N INFERNAL_BRAIN
#PBS -q par16
#PBS -l nodes=1:ppn=16
#PBS -e snakemake_infernal_brain.err
#PBS -o snakemake_infernal_brain.out

source /home/lovelace/proj/proj832/fvperes/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake -p -s Snakefile --cores 16
