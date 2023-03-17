#!/bin/bash

#PBS -N DeepPlnc_GPU
#PBS -q umagpu
#PBS -e DeepPlnc_GPU.err
#PBS -o DeepPlnc_GPU.out

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=16

source /home/lovelace/proj/proj832/fvperes/miniconda3/etc/profile.d/conda.sh
conda activate deepplnc_labis

python3 predict_GPU.py temp
