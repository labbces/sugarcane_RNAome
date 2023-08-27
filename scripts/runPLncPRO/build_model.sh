#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -t 1
#$ -l h=bcmsrv
#$ -tc 1
#$ -pe smp 14

proteincoding=plncpro_data/plant_new_fasta/monocot/train/monocot_pct_train.fa
longnoncoding=plncpro_data/plant_new_fasta/monocot/train/monocot_lnct_train.fa
BLASTDB=/Storage/data1/felipe.peres/swissprot/uniprot/uniprotdb

module load miniconda3
module load blast/2.8.1+
conda activate plncpro

/usr/bin/time -v plncpro build -p $proteincoding -n $longnoncoding -o monocot_model -m monocot.model -d $BLASTDB -t $NSLOTS
