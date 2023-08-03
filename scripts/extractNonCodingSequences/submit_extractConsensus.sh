#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

module load miniconda3
/usr/bin/time -v ./extractConsensusSequences.py -transcriptome 48_transcriptomes.cpc_ncrnas.complete -list putative_ncRNA_consensus.txt

