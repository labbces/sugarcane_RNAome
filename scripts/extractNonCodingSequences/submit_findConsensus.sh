#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -pe smp 1

/usr/bin/time -v ./findConsensus.py CPC2_ncRNAs_list.txt CPC2_PLncPRO_ncRNAs_list.txt RNAplonc_ncRNAs_list.txt putative_ncRNA_consensus.txt
