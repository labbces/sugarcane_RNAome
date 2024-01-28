#!/bin/bash

# -t 1-12
#$ -cwd
#$ -q all.q
#$ -pe smp 1
# -tc 20

module load mcl/14-137

#export m=$((SGE_TASK_ID-1))
#LC_NUMERIC="en_US.UTF-8".
#num=($(LC_NUMERIC="en_US.UTF-8". seq 1.3 0.5 6))

in=Perlo2022_mcl_out.txt
out=${in/.txt}_I2.0.txt

/usr/bin/time -v mcl $in -I 2 -te 1 --abc -o $out
