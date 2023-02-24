#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -t 1-48
#$ -tc 10
#$ -pe smp 1

ASSEMBLY=`ls -1 ../../sugarcaneTranscriptomes/*.fasta | head -n $SGE_TASK_ID | tail -n1`
CPC2=/Storage/data1/felipe.peres/Sugarcane_ncRNA/CPC2_standalone-1.0.1/bin/CPC2.py
ASSEMBLY_BASENAME=`basename $ASSEMBLY`
OUTPUT=${ASSEMBLY_BASENAME/.fasta/_CPC2_output}

module load miniconda3
$CPC2 -i $ASSEMBLY -o ../results/$OUTPUT
