#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -t 1-48
#$ -tc 10
#$ -pe smp 1

input=`ls -1 ../results/*_transcriptome.cpc_ncrnas.fa | head -n $SGE_TASK_ID | tail -n1`
prefixOut=`basename ${input/.fa/_stats}`

module load miniconda3

./getSeqStats.py $input fasta $prefixOut
