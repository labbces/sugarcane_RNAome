#!/bin/bash

#PBS -N MCL
#PBS -q par128
#PBS -l nodes=1:ppn=128
#PBS -m abe
#PBS -e run_mcl_cv1.2.sh.err
#PBS -o run_mcl_cv1.2.sh.out

cd $PBS_O_WORKDIR

source /home/lovelace/proj/proj832/fvperes/miniconda3/etc/profile.d/conda.sh
conda activate MCL

tasks=(1 2 3 4 5 6 7 8)

LC_NUMERIC="en_US.UTF-8"
num=($(LC_NUMERIC="en_US.UTF-8" seq 2.3 0.5 6))

for task in "${tasks[@]}"; do
    export m=$((task-1))

    in=$(ls -1 *_mcl.txt)
    out=pearson_mcl_cv1.2/

    mkdir $out
    mcl $in -I ${num[${m}]} -te 128 -analyze y -o ${out}out.${num[${m}]} --abc

done
