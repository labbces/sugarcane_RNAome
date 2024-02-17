#!/bin/bash

#$ -t 1-11
#$ -cwd
#$ -q all.q
#$ -pe smp 5
#$ -tc 5

module load mcl/14-137

export m=$((SGE_TASK_ID-1))

LC_NUMERIC="en_US.UTF-8"
num=($(LC_NUMERIC="en_US.UTF-8" seq 1.3 0.5 6))

in=`ls -1 *_mcl.txt`
out=pearson_mcl/

mkdir $out
mcl $in -I ${num[${m}]} -te 1 --abc -analyze y -o ${out}out.${num[${m}]} 
