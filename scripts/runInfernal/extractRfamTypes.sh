#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -t 1-48
#$ -tc 1
#$ -pe smp 1

module load miniconda3

cmscan=`ls -1 results/individualGenotypes/*.complete.deoverlapped.cmscan.tblout | head -n $SGE_TASK_ID | tail -n 1`
genotype=`basename ${cmscan/.complete.deoverlapped.cmscan.tblout}`
rfamtypes='rfam-types.txt'
output=results/individualGenotypes/${genotype}.complete.deoverlapped.cmscan.RfamTypes.tblout

/usr/bin/time -v ./extractRfamTypes.py --rfamtypes $rfamtypes --cmscan $cmscan --output $output
