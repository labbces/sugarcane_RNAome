#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load miniconda3

/usr/bin/time -v ./classifyTypeByLength.py -i putative_ncRNAs_length.tsv
