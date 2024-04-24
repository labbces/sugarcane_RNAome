#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -t 1-48
#$ -tc 8
#$ -pe smp 1

# pwd
# /Storage/data1/felipe.peres/Sugarcane/Panzzer/code

GENOTYPE=`ls ../data/ | awk -F "_" '{print $1}' | head -n $SGE_TASK_ID | tail -n1`
GENOTYPE_FILE=${GENOTYPE}_PEP.fix.nr100.fasta
PANZZER=../SANSPANZ.3/runsanspanz.py
SPECIE="Saccharum officinarum x Saccharum spontaneum"

module load miniconda3
module load R/3.5.1

python ${PANZZER} -R -o ",${GENOTYPE}_DE.out,${GENOTYPE}_GO.out,${GENOTYPE}_anno.out" -s "${SPECIE}" < ../data/${GENOTYPE_FILE}
