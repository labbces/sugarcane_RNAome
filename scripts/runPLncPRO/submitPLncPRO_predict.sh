#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -t 1-48
#$ -l h=bcmsrv
#$ -tc 10
#$ -pe smp 1

NONCODING=`ls -1 ../../3_noncodingSequences/results/*.fa | head -n $SGE_TASK_ID | tail -n1`
NONCODING_BASENAME=`basename $NONCODING`

DIAMOND=${NONCODING_BASENAME/.fa/.diamondx}
DIAMOND_OUT=`ls -1 ../../4_DIAMONDx/results/${DIAMOND}`

OUTPUT=${NONCODING_BASENAME/.fa/.plncpro}
MODEL=monocot_model/monocot.model

module load miniconda3
conda activate plncpro

# -r: clean up intermediate files
plncpro predict -i $NONCODING -o ../results/ -p $OUTPUT --blastres $DIAMOND_OUT -t $NSLOTS -m $MODEL -r 
