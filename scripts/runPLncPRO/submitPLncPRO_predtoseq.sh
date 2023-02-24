#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -t 1-48
#$ -l h=bcmsrv
#$ -tc 10
#$ -pe smp 1

NONCODING=`ls -1 ../../3_noncodingSequences/results/*.fa | head -n $SGE_TASK_ID | tail -n1`
NONCODING_BASENAME=`basename $NONCODING`

OUTPUT_PREDICT_FILE=`ls -1 ../results/${NONCODING_BASENAME/.fa/.plncpro}`

OUTPUT_PREDTOSEQ=${OUTPUT_PREDICT_FILE/.plncpro/.lncrnas.plncpro}

module load miniconda3
conda activate plncpro

# -l=0 for lncrnas, -l=1 for mrnas
# -m 200, min_length

#-f ../../3_noncodingSequences/results/B1_transcriptome.cpc_ncrnas.fa
#-o ../results/B1_transcriptome.cpc_ncrnas.lncrnas.plncpro
#-p ../results/B1_transcriptome.cpc_ncrnas.plncpro

plncpro predtoseq -f $NONCODING -o $OUTPUT_PREDTOSEQ -p $OUTPUT_PREDICT_FILE -l 0
