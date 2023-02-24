#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -t 1-48
#$ -tc 10
#$ -pe smp 1

ASSEMBLY=`ls -1 ../../1_sugarcaneTranscriptomes/*.fasta | head -n $SGE_TASK_ID | tail -n1`
CPCFILENAME=`basename ${ASSEMBLY/.fasta/_CPC2_output.txt}`
CPCOUT=`ls -1 ../../2_codingPotential/results/${CPCFILENAME}`

module load miniconda3
/usr/bin/time -v ./extractCPC2NonCodingSequences.py -cpcout $CPCOUT -transcriptome $ASSEMBLY
