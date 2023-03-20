#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -t 1-48
#$ -tc 2
#$ -l h=figsrv
#$ -pe smp 8

module load miniconda3
conda activate DeepPlnc

#Path RNAFold bin
RNAfold=/Storage/data1/felipe.peres/Sugarcane_ncRNA/softwareInstallation/ViennaRNA-2.5.1_gcc7.5.0/bin/RNAfold

#Input transcriptome

FILE_INPUT=`ls -1 ../step1/*fa_mod_single_coding.fa | head -n $SGE_TASK_ID | tail -n1`

BASENAME_INPUT=`basename ${FILE_INPUT}`

#Output from DeepPlnc
TEMP_FILE=${BASENAME_INPUT}_temp

#Bash
$RNAfold --jobs=${NSLOTS} --noPS $FILE_INPUT | paste - - - | awk '{print $3}' > ${BASENAME_INPUT}_RNA_coding

#Bash
cat $FILE_INPUT | paste - - ${BASENAME_INPUT}_RNA_coding > $TEMP_FILE

echo "BASENAME_INPUT =" $BASENAME_INPUT >> jobs2.txt
echo "BASENAME_INPUT_process =" ${FILE_INPUT} >> jobs2.txt
echo "BASENAME_INPUT_RNAfold =" ${BASENAME_INPUT}_RNA_coding >> jobs2.txt
echo "BASENAME_INPUT_final =" $TEMP_FILE >> jobs2.txt
