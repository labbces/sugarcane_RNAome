#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -t 1
#$ -pe smp 1

module load miniconda3
module load salmon/1.8.0
conda activate /home/felipe.peres/.conda/envs/snakemake-tutorial/envs/snakemake_conekt

snakemake -p -k --resources load=10 -s Snakefile --cluster "qsub -q all.q -V -cwd -l h={params.server} -pe smp {threads}" --jobs 10 --latency-wait 60
