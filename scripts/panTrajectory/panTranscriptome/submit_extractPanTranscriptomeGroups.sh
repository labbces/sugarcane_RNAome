#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 1

INFILE=Orthogroups_SugarcanePanTranscriptome_I2.8.tsv

module load Python/3.7.2
/usr/bin/time -v python3.7 computeOrthogroupStats_0.8.py --orthogroupsFile ${INFILE} --numberSpecies 50 --suffixOut I2.8
