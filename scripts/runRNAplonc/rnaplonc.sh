#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load CD-HIT/4.8.1
module load miniconda3
conda activate /home/felipe.peres/.conda/envs/snakemake-tutorial/envs/snakemake_conekt

snakemake -p -s Snakefile --resources load=14 --cluster "qsub -q all.q -V -cwd -pe smp {threads}" --jobs 14 --latency-wait 60
