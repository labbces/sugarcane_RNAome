#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V
#$ -pe smp 1

# create header
echo "clusters,max,ctr,avg,min,DGI,TWI,TWL,sgl,qrt,efficiency,massfrac,areafrac,target-name,source-name" > clminfo.csv

infiles=`ls -1 out*`

# Organizing infos (sep = ',')
for file in $(ls -1 $infiles); do tail -n4 $file | awk '{ for (i=1; i<=NF; i++) { split($i, arr, "="); printf "%s%s%s", sep, arr[2], (i==NF) ? "\n" : "," } }' | tr "\n" ","; echo ""; done >> clminfo.csv
