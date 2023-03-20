#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -t 1-48
#$ -tc 1
#$ -l h=figsrv
#$ -pe smp 1

NONCODING=`ls -1 *.cpc_ncrnas.fa | head -n 41 | tail -n1`

#change header
sed 's/\s.*$//g' $NONCODING > ${NONCODING}_mod

#change sequences
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ${NONCODING}_mod > ${NONCODING}_mod_single

#Loading conda and activating environment
module load miniconda3
conda activate DeepPlnc

# Running DeepPlnc_step1.py
python3 DeepPlnc_step1.py ${NONCODING}_mod_single

# Removing temporary files
rm ${NONCODING}_mod
rm ${NONCODING}_mod_single

# Generating report
echo '### REFORMATING INPUT ###' >> jobs1.txt
echo '# CHANGING HEADER' >> jobs1.txt
echo 'IN: ' $NONCODING >> jobs1.txt
echo 'OUT: ' ${NONCODING}_mod >> jobs1.txt
echo '# CHANGING SEQUENCE FORMAT' >> jobs1.txt
echo 'IN: ' ${NONCODING}_mod >> jobs1.txt
echo 'OUT: ' ${NONCODING}_mod_single >> jobs1.txt
echo '# RUNNING DeepPlnc.py' >> jobs1.txt
echo 'IN: ' ${NONCODING}_mod_single >> jobs1.txt
echo 'OUT: ' ${NONCODING}_mod_single_coding.fa >> jobs1.txt
echo '### COMPLETED ###' >> jobs1.txt
echo ' ' >> jobs1.txt
