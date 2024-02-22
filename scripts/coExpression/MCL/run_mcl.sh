#!/bin/bash

#$ -cwd
#$ -q all.q
#$ -pe smp 10
#$ -t 1-8
#$ -tc 2
#$ -l h=neotera

export m=$((SGE_TASK_ID-1))

LC_NUMERIC="en_US.UTF-8"
# [1.3 1.8 2.3 2.8 3.3 3.8 4.3 4.8 5.3 5.8]
# num=($(LC_NUMERIC="en_US.UTF-8" seq 1.3 0.5 6))

# [1.3 1.8 2.3 2.8 3.3 3.8 4.3 4.8]
num=($(LC_NUMERIC="en_US.UTF-8" seq 1.3 0.5 5))

in=`ls -1 *_mcl.txt`
out=pearson_mcl/

module load mcl/14-137

mkdir $out
/usr/bin/time -v mcl $in -I ${num[${m}]} -te $NSLOTS --abc -analyze y -o ${out}out.${num[${m}]} 
