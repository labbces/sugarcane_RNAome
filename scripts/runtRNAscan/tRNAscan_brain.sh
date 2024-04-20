#!/bin/bash

#PBS -N tRNAscan
#PBS -q par16
#PBS -l nodes=1:ppn=16
#PBS -e tRNAscan_brain.sh.err
#PBS -o tRNAscan_brain.sh.out

cd $PBS_O_WORKDIR

source /home/lovelace/proj/proj832/fvperes/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

snakemake -p -s Snakefile --cores 16 --rerun-incomplete
