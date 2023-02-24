#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -t 1-48
#$ -tc 6
#$ -pe smp 5

UNIREF90=/Storage/data1/felipe.peres/UniRef90/uniref90.dmnd
CPC_NONCODING=`ls -1 ../../3_noncodingSequences/results/*.fa | head -n $SGE_TASK_ID | tail -n 1`
CPC_NONCODING_BASENAME=`basename $CPC_NONCODING`
OUTPUT=${CPC_NONCODING_BASENAME/.fa/.diamondx}

/Storage/data1/felipe.peres/Sugarcane_ncRNA/diamond_2.0.15/diamond blastx --db $UNIREF90 -q ../../3_noncodingSequences/results/$CPC_NONCODING -o ../results/$OUTPUT -p $NSLOTS --outfmt 6 qseqid sseqid pident evalue nident qcovhsp score bitscore qframe qstrand
