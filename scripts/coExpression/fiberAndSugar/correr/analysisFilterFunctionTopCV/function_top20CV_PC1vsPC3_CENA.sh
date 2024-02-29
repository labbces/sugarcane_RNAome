#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -t 1
#$ -tc 1
#$ -pe smp 1

export n=$((SGE_TASK_ID-1))
LC_NUMERIC="en_US.UTF-8"
array=(1 1.3 1.5 1.8 2.0)

module load R/4.0.0

/usr/bin/time -v Rscript filter_VST_function_top20CV_PC1vsPC3.R 2.0
