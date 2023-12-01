#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 5

module load Python/3.7.2

INFILE=`ls -1 *0.txt`
NETWORK=${INFILE/.txt/_network.txt}
MCL=${INFILE/.txt/_mcl_out.txt}

/usr/bin/time -v python3 pcc.py $INFILE $NETWORK $MCL
