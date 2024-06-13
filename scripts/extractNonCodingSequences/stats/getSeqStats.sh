#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

input=putative_ncRNA_consensus.fa
prefixOut=putative

module load miniconda3

./getSeqStats.py $input fasta $prefixOut
