#!/bin/bash

#PBS -N IDENTITY
#PBS -q par128
#PBS -l nodes=1:ppn=128
#PBS -e identity_par128.sh.err
#PBS -o identity_par128.sh.out

cd $PBS_O_WORKDIR

source /home/lovelace/proj/proj832/fvperes/miniconda3/etc/profile.d/conda.sh
conda activate snakemake

#--latency-wait aprox: max secs par128 = 604800. 604800/2304 jobs = 262.5s per job

snakemake -p -s Snakefile --cores 128 --latency-wait 250
